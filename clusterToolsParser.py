from antismash.modules.hmm_detection import rule_parser
from Bio import SearchIO

import string
from enum import IntEnum

class RuleSyntaxError(SyntaxError):
    pass

def is_legal_identifier(identifier):
    """ Returns true if the identifier matches the form:
        [a-zA-Z]{[a-zA-Z0-9_-]}*
    """
    if not identifier[0].isalpha():
        return False
    for char in identifier:
        if not (char.isalpha() or char.isdigit() or char in ['_', '-']):
            return False
    # for now, keep cluster reserved to avoid confusion with previous system
    if identifier == "cluster":
        return False
    # to avoid confusion with minscore
    if identifier == "score":
        return False
    return True

class TokenTypes(IntEnum):
    GROUP_OPEN = 1
    GROUP_CLOSE = 2
    LIST_OPEN = 3
    LIST_CLOSE = 4
    IDENTIFIER = 6
    MINIMUM = 7
    CDS = 8
    AND = 9
    OR = 10
    NOT = 11
    INT = 12
    COMMA = 13
    SCORE = 14

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.name).lower() # the str() is for pylint's sake

    @classmethod
    def classify(cls, text):
        classification = Tokeniser.mapping.get(text)
        if classification is None:
            if text.isdigit() or text[0] == '-' and text[1:].isdigit():
                return cls.INT
            elif is_legal_identifier(text):
                return cls.IDENTIFIER
            else:
                raise RuleSyntaxError("Unclassifiable token: %s" % text)
        return classification

    def is_multichar(self):
        return self.value > 5 # the last single char terminal is 5

class Tokeniser():
    mapping = {"(" : TokenTypes.GROUP_OPEN, ")" : TokenTypes.GROUP_CLOSE,
               "[" : TokenTypes.LIST_OPEN, "]" : TokenTypes.LIST_CLOSE,
               "and" : TokenTypes.AND, "or" : TokenTypes.OR,
               "not" : TokenTypes.NOT, "," : TokenTypes.COMMA,
               "minimum" : TokenTypes.MINIMUM, "cds" : TokenTypes.CDS,
               "minscore" : TokenTypes.SCORE}

    def __init__(self, text):
        """ Processes the given text into a list of tokens """
        self.text = text
        self.tokens = []
        self.tokenise()
        self.current_symbol = []

    def tokenise(self):
        self.current_symbol = []
        position = 0
        while position < len(self.text):
            char = self.text[position]
            # whitespace of any kind separates symbols of interest
            if char in string.whitespace:
                self._finalise(position)
                position += 1
                continue
            # these are their own tokens always, and always single chars
            if char in Tokeniser.mapping:
                self._finalise(position)
                self.tokens.append(Token(char, position))
            # part of a multi-char symbol
            elif char.isalnum() or char in ['-', '_']:
                self.current_symbol.append(char)
            # a comment, don't parse anything in the line
            elif char == "#":
                self._finalise(position)
                try:
                    # skip to after a newline if it exists
                    position += self.text[position:].index('\n') + 1
                except ValueError:
                    # no newline after #, so we're finished
                    return
            else:
                raise RuleSyntaxError("Unexpected character in rule: %s\n%s\n%s^" %(
                        char, self.text, " "*position))
            position += 1
        self._finalise(len(self.text))

    def _finalise(self, position):
        """ convert the current collection of chars into a Token """
        if not self.current_symbol:
            return
        self.tokens.append(Token("".join(self.current_symbol), position))
        self.current_symbol.clear()

class Token():
    def __init__(self, token_text, position):
        self.token_text = token_text
        self.type = TokenTypes.classify(token_text)
        self.position = int(position)
        if self.type.is_multichar():
            self.position -= len(token_text)

    def __getattr__(self, key):
        if key == 'value':
            if self.type != TokenTypes.INT:
                raise AttributeError("Token is not numeric")
            return int(self.token_text)
        elif key == 'identifier':
            if self.type != TokenTypes.IDENTIFIER:
                raise AttributeError("Token has no identifier")
            return self.token_text
        return self.__dict__[key]

    def __repr__(self):
        if self.type == TokenTypes.IDENTIFIER:
            return "'{}'".format(self.token_text)
        if self.type == TokenTypes.INT:
            return self.token_text
        return str(self.type)

class Details():
    def __init__(self, cds, feats, results, cutoff):
        self.cds = cds # str, name of cds that is being classified
        self.features_by_id = feats # { id : feature }
        self.results_by_id = results # { id : HSP list }
        self.possibilities = set(res.query_id for res in results.get(cds, []))
        self.cutoff = cutoff # int

    def in_range(self, cds, other):
        """ returns True if the two Locations are within cutoff distance

            this may be redundant if inputs are already limited, but here
            for safety
        """
        cds_start, cds_end = sorted([cds.start, cds.end])
        other_start, other_end = sorted([other.start, other.end])
        distance = min(abs(cds_end - other_start), abs(other_end - cds_start),
                       abs(cds_end - other_end), abs(other_start - cds_start))
        return distance < self.cutoff

    def just_cds(self, cds_of_interest):
        """ creates a new Details object with cds_of_interest as the focus

            the largest impact is Details.possibilities is updated
        """
        return Details(cds_of_interest, self.features_by_id, self.results_by_id,
                       self.cutoff)

class Conditions():
    def __init__(self, negated, sub_conditions=None):
        self.negated = negated
        assert self.negated in [False, True]
        if sub_conditions is None:
            sub_conditions = []
        self.sub_conditions = sub_conditions
        # just make sure that it's empty or that we have binary ops (a OR b..)
        assert not sub_conditions or len(sub_conditions) % 2 == 1
        if self.sub_conditions:
            operands = self.sub_conditions[::2]
            assert all(isinstance(sub, Conditions) for sub in operands)
            assert all(isinstance(sub, TokenTypes) for sub in self.sub_conditions[1::2])
            for operator in self.sub_conditions[1::2]:
                assert operator in [TokenTypes.AND, TokenTypes.OR]
            unique_operands = set()
            for operand in map(str, operands):
                if operand in unique_operands:
                    raise ValueError("Rule contains repeated condition: %s\nfrom rule %s" %(operand, self))
                unique_operands.add(operand)


    def are_subconditions_satisfied(self, details, any_in_cds=False, local_only=False):
        """
            details : a Details instance
            any_in_cds : bool
            local_only : bool

            any_in_cds is for determining if the CDS being classified contains
                 any hits at all. Effectively changes all conditions to OR, but
                 without looking past CDS being classified

            local_only limits the search to the single CDS in details
        """
        if len(self.sub_conditions) == 1:
            return self.negated ^ self.sub_conditions[0].is_satisfied(details, any_in_cds, local_only)

        # since ANDs are bound together, all subconditions we have here are ORs
        # which means a simple any() will cover us
        sub_results = [sub.is_satisfied(details, any_in_cds, local_only) for sub in self.sub_conditions[::2]]
        return any(sub_results)

    def is_satisfied(self, details, any_in_cds=False, local_only=False):
        return self.negated ^ self.are_subconditions_satisfied(details,any_in_cds, local_only)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        prefix = "not " if self.negated else ""
        if len(self.sub_conditions) == 1 \
                and not isinstance(self.sub_conditions[0], AndCondition):
            return "{}{}".format(prefix, self.sub_conditions[0])
        return "{}({})".format(prefix, " ".join(str(sub) for sub in self.sub_conditions))

class AndCondition(Conditions):
    def __init__(self, subconditions):
        self.operands = subconditions[::2]
        super().__init__(False, subconditions)

    def is_satisfied(self, protein):
        results = [sub.is_satisfied(protein) for sub in self.operands]
        return all(results)

    def __str__(self):
        return " and ".join(map(str, self.operands))

class MinimumCondition(Conditions):
    def __init__(self, negated, count, options):
        self.count = count
        self.options = set(options)
        if len(self.options) != len(options):
            raise ValueError("Minimum conditions cannot have repeated options")
        if count < 1:
            raise ValueError("Minimum conditions must have a required count > 0")
        super().__init__(negated)

    def is_satisfied(self, details, any_in_cds=False, local_only=False):
        hit_count = len(self.options.intersection(set(details.possibilities)))
        if hit_count >= self.count:
            return not self.negated
            # and if we're just checking that one is satisfied, return early too
        if not self.negated and any_in_cds:
            return hit_count > 0
        if local_only:
            return self.negated ^ hit_count >= self.count
        current_cds = details.features_by_id[details.cds]

        # check to see if the remaining hits are in nearby CDSs
        for other_id, other_feature in details.features_by_id.items():
            if other_id == details.cds:
                continue
            if not details.in_range(current_cds.location, other_feature.location):
                continue
            other_options = [r.query_id for r in details.results_by_id.get(other_id, [])]
            hit_count += len(self.options.intersection(set(other_options)))
            if hit_count >= self.count:
                return True ^ self.negated
        return self.negated

    def __str__(self):
        return "{}minimum({}, [{}])".format("not " if self.negated else "",
                self.count, ", ".join(sorted(list(self.options))))

class SingleCondition(Conditions):
    def __init__(self, negated, name):
        self.name = name
        super().__init__(negated)
    def is_satisfied(self, protein):
        # do we only care about this CDS? then use the smaller set

        ## details.possibilities should link to the hmmhit list of the protein
            return self.negated ^ (self.name in protein.getAnnotations('hmm'))
    def __str__(self):
        return "{}{}".format("not " if self.negated else "", self.name)

class ScoreCondition(Conditions):
    def __init__(self, negated, name, score):
        self.name = name
        self.score = score
        super().__init__(negated)

    def is_satisfied(self, details, any_in_cds=False, local_only=False):
        # do we only care about this CDS? then use the smaller set
        if local_only or any_in_cds:
            if self.name in details.possibilities:
                for result in details.results_by_id[details.cds]:
                    if result.query_id == self.name:
                        return self.negated ^ (result.bitscore >= self.score)
            return self.negated

        cds_feature = details.features_by_id[details.cds]
        # look at neighbours in range
        for other, other_hits in details.results_by_id.items():
            other_location = details.features_by_id[other].location
            if not details.in_range(cds_feature.location, other_location):
                continue
            other_possibilities = [res.query_id for res in other_hits]
            if self.name in other_possibilities:
                for result in other_hits:
                    # a positive match, so we can exit early
                    if result.query_id == self.name and result.bitscore >= self.score:
                        return not self.negated

        # if negated and we failed to find anything, that's a good thing
        return self.negated

    def __str__(self):
        return "{}minscore({}, {})".format("not " if self.negated else "", self.name, self.score)

class DetectionRule():
    def __init__(self, conditions):
        self.conditions = conditions

    def detect(self, protein):
        # hit, {hit: 'meh'}, dummyResultsByID, 25
        # change details structure to incorporate protein HMM object
        if not self.conditions.is_satisfied(protein, any_in_cds=True):
            return False #at least one positive match required
        return self.conditions.is_satisfied(protein, any_in_cds=False)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        condition_text = str(self.conditions)
        # strip off outer parens if they exist
        if condition_text[0] == "(" and condition_text[-1] == ')':
            condition_text = condition_text[1:-1]
        return condition_text

class Parser():
    def __init__(self, text):
        self.text = text

        self.current_token = None
        identifiers = set()
        tokens = Tokeniser(text.expandtabs()).tokens
            # gather all signature identifiers
        for token in list(tokens)[1:]:
            if token.type == TokenTypes.IDENTIFIER:
                identifiers.add(token.identifier)
        # start the iterator up for parsing
        self.tokens = iter(tokens)
        try:
            self.current_token = next(self.tokens)
        except StopIteration:
            pass
        self.rule = self._parse_rule()

        #### Add validity check from HMM list. ########
        # verify gathered signature idenfiers exist as signatures
        # unknown = identifiers - signatures.get_signature_names()
        # if unknown:
        #     raise ValueError("Rules contained identifers without signatures: %s" % ", ".join(sorted(list(unknown))))

    def _parse_rule(self):
        """ RULE = classification:ID cutoff:INT extension:INT conditions:CONDITIONS
        """
        conditions = Conditions(False, self._parse_conditions())
        if self.current_token is not None:
            raise RuleSyntaxError("Unexpected symbol %s\n%s\n%s%s" % (
                    self.current_token.type,
                    self.lines[self.current_line - 1],
                    " "*self.current_token.position, "^"))
        return DetectionRule(conditions)

    def _parse_conditions(self):
        """    CONDITIONS = CONDITION {BINARY_OP CONDITIONS}*;
        """

        conditions = []
        lvalue = self._parse_single_condition()
        append_lvalue = True # capture the lvalue if it's the only thing
        while self.current_token and self.current_token.type in [TokenTypes.AND,
                TokenTypes.OR]:
            if self.current_token.type == TokenTypes.AND:
                conditions.append(self._parse_ands(lvalue))
                append_lvalue = False
            else:
                if append_lvalue:
                    conditions.append(lvalue)
                conditions.append(self._consume(TokenTypes.OR))
                lvalue = self._parse_single_condition()
                append_lvalue = True
        if append_lvalue:
            conditions.append(lvalue)

        if self.current_token is None:
            raise RuleSyntaxError("Unexpected end of rule, expected )\n%s" % (
                    self.lines[self.current_line - 1]))
        if self.current_token.type != TokenTypes.GROUP_CLOSE:
            raise RuleSyntaxError("Expected the end of a group, found %s\n%s\n%s^" % (
                    self.current_token,
                    self.lines[self.current_line - 1],
                    " " * self.current_token.position))
        elif self.current_token is not None:
            raise RuleSyntaxError("Unexpected symbol, found %s\n%s\n%s^" % (
                    self.current_token.token_text,
                    self.lines[self.current_line - 1],
                    " " * self.current_token.position))
        return conditions

    def _consume(self, expected):
        if self.current_token is None:
            raise RuleSyntaxError("Unexpected end of line, expected %s" % expected)
        if self.current_token.type != expected:
            raise RuleSyntaxError("Expected %s but found %s\n%s\n%s%s" % (
                    expected, self.current_token.type,
                    self.lines[self.current_line - 1],
                    " "*self.current_token.position, "^"))
        if expected == TokenTypes.INT:
            ret = self.current_token.value
        elif expected == TokenTypes.IDENTIFIER:
            ret = self.current_token.identifier
        else:
            ret = self.current_token.type
            assert isinstance(ret, TokenTypes)
        try:
            self.current_token = next(self.tokens)
        except StopIteration:
            self.current_token = None
        return ret

    def _is_not(self):
        """ [UNARY_OP]
        """
        negated = self.current_token.type == TokenTypes.NOT
        if negated:
            self._consume(TokenTypes.NOT)
        return negated

    def _parse_ands(self, lvalue):
        """ CONDITION and CONDITION { and CONDITION}
            ^ lvalue being passed in
        """
        and_conditions = [lvalue]
        and_conditions.append(self._consume(TokenTypes.AND))
        and_conditions.append(self._parse_single_condition())
        while self.current_token and self.current_token.type == TokenTypes.AND:
            and_conditions.append(self._consume(TokenTypes.AND))
            next_condition = self._parse_single_condition()
            and_conditions.append(next_condition)
        return AndCondition(and_conditions)

    def _parse_single_condition(self, allow_cds):
        """
            CONDITION = [UNARY_OP] ( ID | CONDITION_GROUP | MINIMUM | CDS );
            or we're in a CDS (i.e. allow_cds == False)
            CDS_CONDITION = [UNARY_OP] ID {BINARY_OP CDS_CONDITION}*;

        """
        negated = self._is_not()
        if self.current_token is None:
            raise RuleSyntaxError("Rules cannot end in not")
        if self.current_token.type == TokenTypes.GROUP_OPEN:
            return Conditions(negated, self._parse_group())
        elif self.current_token.type == TokenTypes.MINIMUM:
            return self._parse_minimum(negated=negated)
        elif self.current_token.type == TokenTypes.SCORE:
            return self._parse_score(negated=negated)
        return SingleCondition(negated, self._consume(TokenTypes.IDENTIFIER))

    def _parse_score(self, negated=False):
        """
            SCORE = minscore GROUP_OPEN ID COMMA INT GROUP_CLOSE
            e.g. minscore(trsC, 150)
        """
        self._consume(TokenTypes.SCORE)
        self._consume(TokenTypes.GROUP_OPEN)
        ident = self._consume(TokenTypes.IDENTIFIER)
        self._consume(TokenTypes.COMMA)
        score = self._consume(TokenTypes.INT)
        self._consume(TokenTypes.GROUP_CLOSE)
        return ScoreCondition(negated, ident, score)

    def _parse_group(self, allow_cds):
        """
            CONDITION_GROUP = GROUP_OPEN CONDITIONS GROUP_CLOSE;
        """
        self._consume(TokenTypes.GROUP_OPEN)
        conditions = self._parse_conditions()
        self._consume(TokenTypes.GROUP_CLOSE)
        return conditions

    def _parse_minimum(self, negated=False):
        """
            MINIMUM = MINIMUM_LABEL GROUP_OPEN
                  count:INT COMMA
                  LIST COMMA
                  GROUP_CLOSE;
        """
        initial_token = self.current_token
        self._consume(TokenTypes.MINIMUM)
        self._consume(TokenTypes.GROUP_OPEN)
        count = self._consume(TokenTypes.INT)
        self._consume(TokenTypes.COMMA)
        options = self._parse_list()
        self._consume(TokenTypes.GROUP_CLOSE)
        if count < 0:
            raise ValueError("minimum count must be greater than zero: \n%s\n%s^" % (
                             self.lines[self.current_line - 1],
                             " "*initial_token.position))
        return MinimumCondition(negated, count, options)

    def _parse_list(self):
        """
            LIST = LIST_OPEN contents:COMMA_SEPARATED_IDS LIST_CLOSE;
        """
        self._consume(TokenTypes.LIST_OPEN)
        contents = [self._consume(TokenTypes.IDENTIFIER)]
        while self.current_token.type == TokenTypes.COMMA:
            self._consume(TokenTypes.COMMA)
            contents.append(self._consume(TokenTypes.IDENTIFIER))
        self._consume(TokenTypes.LIST_CLOSE)
        return contents


class hmmHit():
    def __init__(self,hitID):
        self.query_id = hitID

class Details():
    def __init__(self, cds, feats, results, cutoff):
        self.cds = cds # str, name of cds that is being classified
        self.features_by_id = feats # { id : feature }
        self.results_by_id = results # { id : HSP list }
        self.possibilities = set(res.query_id for res in results.get(cds, []))
        self.cutoff = cutoff # int

    def in_range(self, cds, other):
        """ returns True if the two Locations are within cutoff distance

            this may be redundant if inputs are already limited, but here
            for safety
        """
        cds_start, cds_end = sorted([cds.start, cds.end])
        other_start, other_end = sorted([other.start, other.end])
        distance = min(abs(cds_end - other_start), abs(other_end - cds_start),
                       abs(cds_end - other_end), abs(other_start - cds_start))
        return distance < self.cutoff

    def just_cds(self, cds_of_interest):
        """ creates a new Details object with cds_of_interest as the focus

            the largest impact is Details.possibilities is updated
        """
        return Details(cds_of_interest, self.features_by_id, self.results_by_id,
                       self.cutoff)

def prepareRulesByIdDict(hmmResults,cutoff):
    results_by_id = dict()
    for result in hmmResults:
        for hsp in result.hsps:
            if hsp.bitscore > cutoff:
                if hsp.hit_id not in results_by_id:
                    results_by_id[hsp.hit_id] = [hmmHit(hsp.query_id)]
                else:
                    results_by_id[hsp.hit_id].append(hmmHit(hsp.query_id))
    return results_by_id


test = Tokeniser('transatpks	20	20	cds(minscore(PKS_KS,20) and ATd and (PKS_KS or ene_KS or mod_KS or hyb_KS or itr_KS or tra_KS))')
test_iter = iter(test.tokens)
current = next(test_iter)
print(current)
print(test.tokens)
# parser = rule_parser.Parser(open('/Users/emzodls/antismash5/antismash/modules/hmm_detection/test_rule.txt'))
# results = SearchIO.parse("/Volumes/Data/clusterToolsDB/mibig/mibigNRPS.out",'hmmsearch3-domtab')
# dummyResultsByID = prepareRulesByIdDict(results,25)
#
# for hit in dummyResultsByID:
#     details = Details(hit,{hit:'meh'}, dummyResultsByID, 25)
#     if parser.rules[0].conditions.is_satisfied(details):
#         print(hit)
