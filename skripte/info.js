$(function()
{
	$("#infoButton")
       .text("")
       .css(
	{ "z-index":"2",
	  "background":"rgba(0,0,0,0)", "opacity":"0.9",
	  "position":"absolute", "top":"4px", "left":"4px"
	})
       .append("<img width='32' height='32' src='skripte/info.png'/>")
       .button()
	.click(
		function()
		{
			$("#infoBox").dialog("open");
		});

        $("#infoBox")
	.css(
	{
	   "background":"rgba(255,255,255,0.5)"
	})

	.dialog({ autoOpen: true,
		show: { effect: 'fade', duration: 500 },
		hide: { effect: 'fade', duration: 500 },
		draggable: false,
		resizable: false,
		position: {my: "left top", at: "right bottom", of: "#infoButton"},
		open: function() { $(this).parents('.ui-dialog').attr('tabindex', -1)[0].focus(); }
	});
});
