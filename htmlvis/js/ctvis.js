
function ct(ct_data, ct_scale, hits,hit_colors,hmms,hmm_colors,container) {
	var ct_data = ct_data; 
    var hits = hits;
    var hit_colors = hit_colors;
    var hmms = hmms;
    var hmm_colors = hmm_colors;
    
  $("#" + container).html("");
  var ui = $("#" + container);
    ui.append("<h3>ClusterTools Hits</h3>")
  for (var i in ct_data) {
  	ui.append("<h4>" + ct_data[i]["id"] + " (Location: " + ct_data[i]["start"] + "-" + ct_data[i]["end"] + ", Size: " + ct_data[i]["size"] + " nts)" + "<h4>");
    var svg = $(Arrower.drawClusterToolsSVG(ct_data[i],ct_scale));
    svg.css("clear", "both");
    svg.addClass("arrower-svg");
    ui.append(svg);
    ui.find("svg.arrower-svg").each(
    	function(){
    		$(this).attr("width",$(this).find("g")[0].getBoundingClientRect().width);
    		$(this).attr("height",$(this).find("g")[0].getBoundingClientRect().height*1.3);
    	});
  }
    ui.append("<h3>Legend</h3>");
    ui.append("<h4>Queries<h4>");
    ui.append("<div>")
    for (var i in hits){
        ui.append('<svg width="20px" height="20px"><rect width="20px" height="20px" style="fill:' + hit_colors[i]+ ';stroke-width:1;stroke:rgb(0,0,0)"/></svg>')
        ui.append('\t'+hits[i]+'\t')
    }
    ui.append("</div>")
    ui.append("<h4>HMMs<h4>");
    ui.append("<div>")
    for (var i in hmms){
        ui.append('<svg width="20px" height="20px"><rect width="20px" height="20px" style="fill:' + hmm_colors[i]+ ';stroke-width:1;stroke:rgb(0,0,0)"/></svg>')
        ui.append('\t'+hmms[i]+'\t')
    }
    ui.append("</div>")
}