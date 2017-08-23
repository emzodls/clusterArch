
function ct(ct_data, ct_scale, hits,hit_colors,hmms,hmm_colors,container) {
	var ct_data = ct_data; 
    var hits = hits;
    var hit_colors = hit_colors;
    var hmms = hmms;
    var hmm_colors = hmm_colors;
    
  $("#" + container).html("");
  var ui = $("#" + container);
    ui.append("<h3>ClusterTools Hits</h3>");
    ui.append('<div class="row"><div class="col-xl" id="clusters">');
  var clusterContainer = $("#clusters");
  
  for (var i in ct_data) {
    clusterContainer.append('<div class="row"><div class="col-sm-auto">')
    var svgContainer = ui.find(".col-sm-auto").last();
  	svgContainer.append("<h4>" + ct_data[i]["id"] + " (Location: " + ct_data[i]["start"] + "-" + ct_data[i]["end"] + ", Size: " + ct_data[i]["size"] + " nts)" + "</h4>");
    var svg = $(Arrower.drawClusterToolsSVG(svgContainer,ct_data[i],ct_scale));
    svg.css("clear", "both");
    svg.addClass("arrower-svg");
    svgContainer.append(svg);
    svgContainer.find("svg.arrower-svg").each(
    	function(){
    		$(this).attr("width",$(this).find("g")[0].getBoundingClientRect().width);
    		$(this).attr("height",$(this).find("g")[0].getBoundingClientRect().height*1.3);
    	});
  }
    clusterContainer.parent().append('<div class="col sticky-top" id="legend">');
    var legendDiv = $("#legend");
    legendDiv.append('<button class="btn-info pull-right" data-toggle="collapse" data-target="#content">Legend</button>')
    legendDiv.append('<div id="content" class="collapse">')
    var legendContent = $("#content")
    legendContent.append('<div class="row"><h4>Queries<h4></div>');
    legendContent.append('<div class="row">');
    for (var i in hits){
        var legendContainer = ui.find(".row").last();
        legendContainer.append('<div class="col">');
        legendContainer.append('<svg width="30px" height="20px"><rect width="20px" height="20px" style="fill:' + hit_colors[i]+ ';stroke-width:1;stroke:rgb(0,0,0)"/></svg>');
        legendContainer.append('\t'+hits[i]+'\t<br>');
    }
    
    legendContent.append('<div class="row"><h4>HMMs<h4></div>');
    legendContent.append('<div class="row">');
    for (var i in hmms){
        var legendContainer = ui.find(".row").last();
        legendContainer.append('<div class="col">');
        legendContainer.append('<svg width="30px" height="20px"><rect width="20px" height="20px" style="fill:' + hmm_colors[i]+ ';stroke-width:1;stroke:rgb(0,0,0)"/></svg>');
        legendContainer.append('\t'+hmms[i]+'\t<br>');
    }
    ui.append("</div>")
}