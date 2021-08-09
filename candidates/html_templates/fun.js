$(document).ready(function() {
	$( ".candidate" ).hide();
	var fragment_link = window.location.hash.substr(1);
	if(fragment_link != ""){
		$("#" + fragment_link).show(); 
		$("#" + fragment_link).children("figure").each(function(){
			source = $(this).children("img").attr("secret_src");
			$(this).children("img").attr("src", source);
		}).get;
	}
	$( "cand" ).click(function (e) {
		$( ".candidate" ).hide();
		$("#content_" + e.target.id).show(); 
		$("#content_" + e.target.id).children("figure").each(function(){
			source = $(this).children("img").attr("secret_src");
			$(this).children("img").attr("src", source);
		}).get;
	});
});
