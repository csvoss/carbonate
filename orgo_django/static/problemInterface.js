

var jq = jQuery.noConflict();

var isSuccess = false;

var moleculeVisible = false;

function setup(typeableReagents) {
    //For making molecules and reactions draggable
    //$( ".molecule" ).draggable({helper: "clone", revert:true, revertDuration: 100});
    //$( ".reaction" ).draggable({helper: "clone", revert:true, revertDuration: 100});
    $( ".reagent" ).draggable({helper: "clone", revert:true, revertDuration: 100});
    
    isSuccess = false;
    moleculeVisible = false;
    
    $(".helperPopup").draggable();
    
    $("#leftbar").droppable({
        accept: '.reagent',
        drop: function(event, ui) {
            if (tutorial == 2) {
                tutorial = 3;
                $(".wave2").css("display", "none");
                $(".wave3").css("display", "block");
            }
            $.ajax({
                type: "POST",
                url: "/orgo/api/checkSingleStepReaction/",
                data: {'reagents': ui.draggable.attr("reagentString")},
                success: function(data) {
                
                    dataObject = jQuery.parseJSON(data);
                    //Update with new reaction step
                    $("#reactionsHere").html("");
                    $("#reactionsHere").append(ui.draggable.html());
                    //Update with new product
                    //dataObject.product should be svg data
                    $("#productMolecule").html(dataObject.product);
                    
                    
                    //Update with whether or not the user was successful
                    if (dataObject.success) {
                        isSuccess = true;
                        $("#successbox").html("<div style=\"color:#00FF00\" class=\"headnavbar\">SUCCESS!</div>");
                        $("#reagentsHere").html("<div style=\"border-radius: 10px; padding:20px; background-color: #181818; margin-top:50px; font-color:#00FF00\"><h3 style='color:#999999'>Congrats!</h3> <br /> <a href='/orgo/namereagent/'> Do another </a><br /><a href='/orgo/'> Back to home </a></div>");
                        //$("#messageArea").html("<h2 style='color:green'>Congrats!</h2> <br /> <a href='/orgo/namereagent/'> Do another </a><br /><a href='/orgo/'> Back to home </a>");
                        //$("#messageBox").css('display', "block");
                        
                    } else {
                        $("#successbox").html("<div style=\"color:#FFFFFF\" class=\"headnavbar\">Try again.</div>");
                    }
                },
            });
        }
    });
    
    
   
    //For having an autocomplete box which can take in multiple values
    function split( val ) {
        return val.split( /,\s*/ );
    }
    function extractLast( term ) {
        return split( term ).pop();
    }
    //Autocomplete box for reagents. Borrowed from jquery.
    $( "#reagentTyperBox" ).bind( "keydown", function( event ) {
        if ( event.keyCode === $.ui.keyCode.TAB && $( this ).data( "autocomplete" ).menu.active ) {
            event.preventDefault();
        }
    })
    .autocomplete({
        minLength: 0,
        source: function( request, response ) {
            response( $.ui.autocomplete.filter(typeableReagents, extractLast( request.term ) ) );
        },
        focus: function() {
            return false;
        },
        select: function( event, ui ) {
            var terms = split( this.value );
            terms.pop();
            terms.push( ui.item.value );
            terms.push( "" );
            this.value = terms.join( ", " );
            return false;
        }
    });
      
    
    
    //For reagents in the sidebar
    $( "#reagentsHere" ).sortable();
    $( "#reagentsHere" ).disableSelection();

    
    
    
    updateReagents = function() {
        //Don't need to send anything back!
        //Update the reagents present in the sidebar with a new reagent.
        //Keep data about that reagent in currentReagents.
        if (tutorial == 1){
            tutorial = 2;
            $(".wave1").css("display", "none");
            $(".wave2").css("display", "block");
        }
        var reagentString = $("#reagentTyperBox").val();
        
        if (reagentString != "") {
            $.ajax({
                type: "POST",
                url: "/orgo/api/returnReagentHtml/",
                data: {'reagentString': reagentString},
                success: function(data) {
                    if (data != "") {
                        if (isSuccess)
                            $("#reagentsHere").html("")
                        $("#reagentsHere").prepend(data);
                    }
                },
            });
        }
    }

    
    //Make big molecules.
    $("#bigMolecule").each(function(){
        $(this).click(function(){
            $(this).css("left", "-9999px");
            moleculeVisible = false;
        });
    });
    
    updateBigMolecule();
    
    //Submit reagent when enter is pressed
    $("#reagentTyperBox").keydown(function(e){
        if (e.keyCode == 13) {
            updateReagents();
            $("#reagentTyperBox").val("");
            $("#addReagent").focus();
        }
    });
}

//Makes the bigMolecule svg show up when you hover over the corresponding regular molecule.
function updateBigMolecule(){
    $(".molecule, #target").each(function(){
        var bigSelector = "#bigMolecule";
        $(this).click(function(){
            //What happens when mouse enters area
            //Wait for a few seconds, then show the big molecule
            if ($(this).html().length < 6){
                //Sorta hackish way of testing whether anything is in the molecule div
                return
            }
            out=$(this).html().replace('height="200px"', 'height="400px"').replace('width="200px"', 'width="400px"');
            image = "<img src=\"/orgo/static/x.png\" style=\"position:absolute;\"></img>";
            $(bigSelector).css('left', '400px');
            out = image + out;
            $(bigSelector).html(out);
            moleculeVisible = true;

        });
    });
}
    
showAnswer = function(){
    //If user clicks the answer button, make a AJAX request and get the answer.
    $.ajax({
        type: "GET",
        url: "/orgo/api/showSingleStepAnswer/",
        success: function(data) {
            $("#reactionsHere").html(data);
            $("#productMolecule").html($("#target").html());
        }
    });
}

function closeTutorial() {
    tutorial = 0;
    $(".wave3").css("display", "none");
}

var updateReagents;



