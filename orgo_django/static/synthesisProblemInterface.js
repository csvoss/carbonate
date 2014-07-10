
var jq = jQuery.noConflict();
var state = "Normal";
var pk = 0; //Gets changed before chat starts.
var REFRESH = 1500; //Chat refresh time, in ms.
var isSolution = false;
var isSuccess = false;

//For drawing stuff upon initialization
var redrawProblem = function() {
    $.ajax({
        type: "GET",
        url: "/orgo/api/getSynthesisData/",
        data: {},
        success: function(data) {
            drawAllTheThings(data);
        },
    });
}

var successUpdate = function(success) {
    if (success == "solution"){
        $("#successbox").html("<div style=\"color:#EEEEEE\" class=\"headnavbar\">Solution</div>");
    }else if (success){
        onSuccess();
        isSuccess = true;
    }else{
        $("#successbox").html("<div style=\"color:#EEEEEE\" class=\"headnavbar\">Unsolved</div>");
    }
}

var drawAllTheThings = function(data) {
    dataObject = jQuery.parseJSON(data);
    //Redraw all the things
    
    
    $("#leftbar").html("");
    
    try {
        drawMolecules(moleculeListSort(dataObject.molecules, dataObject.arrows));
    } catch(e) {
        //console.log("ERROR! In drawMolecules");
    }
    try {
        jsPlumb.reset();
    } catch(e) {
        ////console.log("ERROR! In jsPlumb reset");
    }
    try {
        drawArrows(dataObject.arrows);
    } catch(e) {
        ////console.log("ERROR! In drawArrows");
    }
    try {
        jsPlumb.repaintEverything();
    } catch(e) {
        ////console.log("ERROR! In jsPlumb.repaintEverything()" + e.stack);
    }
    
    //Update with whether or not the user was successful
    successUpdate(dataObject.success);
}
    
//moleculesSorted is a "sorted" array of arrays, which optimizes arrow flow:
//        [ [[id,svg], [id,svg], [id,svg], ... ], //first row
//        [[id,svg], [id,svg], [id,svg], ... ], ... ] //second row
var drawMolecules = function(moleculesSorted) {
    ////console.log(moleculesSorted.length);
    //we need to create a bunch of divs of the form
        //<div id={{id_number}} class="molecule">{{svg_data}}</div>
    htmlToAddToChart = "";
    for (var i=0; i<moleculesSorted.length; i++){      //for row in moleculesSorted
        htmlToAddToChart += "<div id=\"cleared\">";
        for (var j=0; j<moleculesSorted[i].length; j++) {      //for molecule in row
            molecule = moleculesSorted[i][j];
            var classString = "' class='molecule'>";
            if (molecule[3]) {
                classString = "' class='molecule starting'>";
            }  
            //add a new div
            div = "<div id='" +
                    String(molecule[0]) + 
                    classString +
                    String(molecule[1]) + 
                    "</div>";
            htmlToAddToChart += div;
            ////console.log("Added new div.");
        }
        //add a line break, or something
        htmlToAddToChart += "</div><br />";
    }
    
    
    $(".molecule").remove();
    //put the constructed html in #leftbar
    $("#leftbar").html("");
    $("#leftbar").append(htmlToAddToChart);
    //For making molecules and reactions draggable
    $( ".molecule" ).draggable({helper: "clone", revert:"invalid", revertDuration: 100, drop: function(event, ui) {} });
    //The right bar is droppable, and triggers deletion of molecules
    $("#rightbar, #offscreen").droppable({
        greedy: true,
        tolerance: 'pointer',
        drop: function(event, ui) {
            
            //draw things
            if (ui.draggable.hasClass("molecule") && !(ui.draggable.hasClass("starting")) && !isSolution) {
                try {
                    id = ui.draggable.attr("id");
                    $.ajax({
                        type: "POST",
                        url: "/orgo/api/deleteMolecule/",
                        data: {"moleculeID": id},
                        success: function(data) {
                            ui.draggable.remove();
                            drawAllTheThings(data);
                            ////console.log("Post request successful.");
                        }
                    });
                } catch(e) {
                    ////console.log("ERROR! In rightbar/offscreen droppable.");
                }
            }
        }
    });
    $("#leftbar").droppable({
        greedy: true,
        drop: function(event, ui) {
        }
    });
    //Molecules are droppable, triggering addition of reagents or other molecules
    $(".molecule").droppable({
        greedy: true,
        accept: '.molecule, .reagent',
        hoverClass: 'drophover',
        drop: function(event, ui) {
            if (isSolution) {
            /*
                isSolution = false;
                redrawProblem();
                try{
                    jsPlumb.repaintEverything();
                } catch(e) {
                    ////console.log("ERROR! 42 In jsPlumb.repaintEverything();");
                }
                */
            }
            else if (ui.draggable.hasClass("molecule")) {
                $.ajax({
                    type: "POST",
                    url: "/orgo/api/addMoleculeToMolecule/",
                    data: {'molecule1': ui.draggable.attr("id"),
                           'molecule2': $(this).attr("id")},
                    success: function(data) {
                        drawAllTheThings(data);
                    },
                });
            }
            else if (ui.draggable.hasClass("reagent")) {
                if (tutorial == 2) {
                    tutorial = 3;
                    $(".wave2").css("display", "none");
                    $(".wave3").css("display", "block");
                }
                $.ajax({
                    type: "POST",
                    url: "/orgo/api/addReagentToMolecule/",
                    data: {'reagents': ui.draggable.attr("reagentString"),
                           'moleculeOn': $(this).attr("id")},
                    success: function(data) {
                        drawAllTheThings(data);
                        ////console.log("Reagent string sent to ajax request: "+ui.draggable.attr("reagentString"));
                    },
                });
            }
        }
    });
    updateBigMolecule();
}


//molecules is an array of arrays: [ [idnumber, "<svg>...</svg>"], ... ]
//arrows is an array of arrays: [ [idnumber1, idnumber2, "reagentText"], ...]
//returns a "sorted" array of arrays, which optimizes arrow flow:
//        [ [[id,svg], [id,svg], [id,svg], ... ], //first row
//        [[id,svg], [id,svg], [id,svg], ... ], ... ] //second row
//CAN BE MADE BETTER, specifically by minimizing arrows crossing over each other.
//RETURN TO THIS LATER.
var moleculeListSort = function(molecules, arrows) {
    ////console.log("Init molecule len: "+String(molecules.length));
    //Iterate through the list, determining horizontal rank:
        //Figure out which molecules are first (aka starting molecules)
    
    // This step is O(n^2).
    // It *could* be O(n).
    
    var arrowProducts = [];
    for (var i=0; i<arrows.length; i++) {
        arrowProducts.push(arrows[i][1]);
        ////console.log("Arrow: "+arrows[i][0]+", "+arrows[i][1]+", "+arrows[i][2]);
    }
    
    
    var startingMolecules = [];
    for (var i=0; i<molecules.length; i++) {
        m = molecules[i];
        if (arrowProducts.indexOf(m[0]) == -1) {
            startingMolecules.push([m[0], m[1], 0, true]); //the "true" is to keep track of it being a starting molecule.
            //console.log("Added to starting: "+m[0]);
        }
        //else console.log("Not added to starting: "+m[0]);
    }
    
    
    //Figure out the minimum distance of all other molecules from the starting molecules
    //Recursion!
    var svgGet = function(idNumber) {
        for (var i = 0; i<molecules.length; i++)
            if (molecules[i][0] == idNumber)
                return molecules[i][1];
        return "null";
    }
    
    //This entire next block is devoted to finding out the row (an integer) at which to draw a molecule.
    var maxInd = 0;
    var rank = function(currentMolecules, ind) {
        //If the size of currentMolecules is equal to the size of molecules, return -- you're done
        if (currentMolecules.length >= molecules.length)
            return currentMolecules;
            
            
        maxInd = ind;
        //Go through all of the reactions with beginnings in currentMolecules.
        //If not already in currentMolecules, append them, with the relevant index of distance from starting, to a new list
        //Call rank() again on the new list, with an incremented index.
        //This is O(n^3), and it also has to be called multiple times via recursion >.<            
        
        //[x[0] for each (x in currentMolecules)]
        
        var currentMoleculeIndices = []
        //var s = ""
        for (var j=0; j<currentMolecules.length; j++) {
            currentMoleculeIndices.push(currentMolecules[j][0]);
        //    s +=  currentMolecules[j][0]+",";
        }
        
        //[ [arrow[0], svgGet(arrow[0]), ind] for each (arrow in arrows) if ([x[0] for each (x in currentMolecules)].indexOf(arrow[0])==-1)  ]
        var toAdd = []
        for (var i=0; i<arrows.length; i++) {
            var aP = arrows[i][1]; //arrow product
            var aR = arrows[i][0]; //arrow reactant
            
            var toAddContents = []
            for (var j=0; j<toAdd.length; j++) {
                toAddContents.push(toAdd[j][0]);
            }

            
            //If the current arrow's product is not in the current list AND not in toAdd[any][0],
            //AND the current arrow's reactant IS in the current list,
            //add it
            if ((currentMoleculeIndices.indexOf(aP) == -1) && (toAddContents.indexOf(aP) == -1) && !(currentMoleculeIndices.indexOf(aR) == -1))
                toAdd.push([aP, svgGet(aP), ind, false]);
            //else
            //    console.log("Arrow product "+aP+" is not in current list, "+s);
        }
        
        if (ind > 10)
            return currentMolecules;
        
        //console.log("Re-recursing at 1+"+ind);
        return rank(currentMolecules.concat(toAdd), ind+1);
    }
    
    //console.log("Starting recursion.");
    //rankedMolecules is an array of arrays: [ [idnumber, "<svg>...</svg>", distanceFromStartingMolecules], ... ]
    var rankedMolecules = rank(startingMolecules, 1);
    //console.log("Ended recursion.");
    
    //console.log("Intermediate molecule len: "+String(rankedMolecules.length));
    s = "";
    for (var i = 0; i<rankedMolecules.length; i++)
        s += rankedMolecules[i][0] + ", " + rankedMolecules[i][2] + "\n";
    //console.log(s); 
       
    //Iterate through the list, assigning nodes to vertical columns
    //Produce the final list.
    output = [startingMolecules];
    for (var i = 1; i <= maxInd; i++) {
    
        //[ [ triple for each (triple in rankedMolecules) if (triple[2] == i) ] ]
        //alsoToAdd represents the next row to add to the output array.
        alsoToAdd = [];
        for(var j = 0; j < rankedMolecules.length; j++) {
            if (rankedMolecules[j][2] == i) {
                //console.log("...");
                //console.log(rankedMolecules[j]);
                alsoToAdd.push(rankedMolecules[j]);
                //console.log("!!!");
                //console.log(alsoToAdd);
            }
        }
        
        //This bit of code: Rearranges the items in alsoToAdd so as to minimize arrow-crossing.
        var prevToAdd = output[output.length-1]; //the previous row
        var finalToAdd = [];
        for (var j = 0; j < prevToAdd.length; j++) {
            //Add, to finalToAdd, (the kth) item in alsoToAdd if it is at the end of any arrow starting with (the jth) element in prevToAdd
            //This is kind of like a sort, so it might be time consuming. But we're O(n^3) from earlier, anyways. >.>
            
            var arrowProdsToCheck = []; //array will contain all products of arrows beginning with the jth element in prevToAdd
            for (var l = 0; l < arrows.length; l++) {
                if (arrows[l][0] == prevToAdd[j][0]) { //if the reactant is in the previous row, allow the product in array
                    arrowProdsToCheck.push(arrows[l][1]);
                }
            } //end initialization of arrowProdsToCheck
            
            for (var k = 0; k < alsoToAdd.length; k++) {
                //Check if: the kth item in alsoToAdd is --
                    //not already in finalToAdd (automatic no) --this is important! some molecules are products of two arrows.
                    //at the end of any arrow in arrowProdsToCheck (automatic yes, UNLESS the above)
                    //if not any of the above, then default to don't add it
                var shouldAddKthItem = false;
                if (arrowProdsToCheck.indexOf(alsoToAdd[k][0]) != -1) { //it's in arrowProdsToCheck, so probably add it
                    shouldAddKthItem = true;
                    if (finalToAdd.indexOf(alsoToAdd[k]) != -1) { //check if it's already in finalToAdd, though; if it is, DON'T add it
                        shouldAddKthItem = false;
                    }
                } 
                
                if (shouldAddKthItem) {
                    finalToAdd.push(alsoToAdd[k]);
                }
            } 
        }
        //Add, to finalToAdd, anything in alsoToAdd which is NOT yet in finalToAdd
        for (var k = 0; k < alsoToAdd.length; k++) {
            if (finalToAdd.indexOf(alsoToAdd[k]) == -1) {
                finalToAdd.push(alsoToAdd[k]);
            }
        }
        
        //output.push(alsoToAdd); //in case of bugs with finalToAdd, uncomment this
        output.push(finalToAdd);
    }
    
    //console.log("Finit molecule height: "+String(output.length));
    //for(var i = 0; i<output.length; i++)
        //console.log("Row length: "+String(output[i].length));
    
    
    
    
    
    
    
    return output
}



//arrows is a list of lists: [ [idnumber1, idnumber2, "reagentText"], ...]
var drawArrows = function(arrows) {
    //clear existing jsplumb connections -- how to?
    
    jsPlumb.bind("ready", function() {
        try {
        jsPlumb.setSuspendDrawing(true);
        } catch(e) {
        //console.log("ERROR 52");
        }
        for (var i = 0; i < arrows.length; i++) {
            try {
            molecule1 = document.getElementById(String(arrows[i][0]));
            molecule2 = document.getElementById(String(arrows[i][1]));
            } catch(e) {
            //console.log("ERROR 53");
            }
            try {
            var conn = jsPlumb.connect({
                source:molecule1,  // just pass in the current node in the selector for source 
                target:molecule2,
                // here we supply a different anchor for source and for target, and we get the element's "data-shape"
                // attribute to tell us what shape we should use, as well as, optionally, a rotation value.
                anchors:[
                    [ "Perimeter", {shape:"rectangle"}],
                    [ "Perimeter", {shape:"rectangle"}]
                ],
                overlays:[
                    [ "Label", {label:((arrows[i][2] == "")? (""):("<span class=\"arrowLabel\">"+arrows[i][2]+"</span>")), id:"label"}]
                ],
                enabled:false,
                endpoint:["Dot", {radius:1}],
            });
            } catch(e) {
            //console.log("ERROR 54");
            }
            ////console.log("...");
            //conn.getOverlay("label").setLabel(conn.getParameters().reagents);
        }
        try {
        jsPlumb.setSuspendDrawing(false);
        } catch(e) {
        //console.log("ERROR 55");
        }
    });
}





//For having an autocomplete box which can take in multiple values
function split( val ) {
    return val.split( /,\s*/ );
}
function extractLast( term ) {
    return split( term ).pop();
}


//Start-up (document.ready) functions run by both the solver and the helper.
function allSetup() {
    try {
        jsPlumb.bind("ready", function() {
        
            $(".helperPopup").draggable();
            
            $('#leftbar').scroll(function () {
                try {
                    jsPlumb.repaintEverything();
                } catch(e) {
                    //console.log("ERROR! 43 In jsPlumb.repaintEverything");
                }
            });

            try {
                jsPlumb.importDefaults({
                    Connector:"StateMachine",
                    PaintStyle:{ lineWidth:3, strokeStyle:"#181818"},
                    ConnectionOverlays : [
                        [ "Arrow", { 
                            location:1,
                            id:"arrow",
                            length:10,
                            foldback:1
                        } ],
                    ],
                    Anchors : [ "TopCenter", "BottomCenter" ]
                });
            } catch(e) {
                //console.log("ERROR! 43 In jsPlumb.importDefaults;");
            }
        });
    } catch(e) {
        //console.log("ERROR! In allSetup() ???");
    }
    
    //On resize, recalculate the two div heights.
    $(window).resize(function(){
        remakeDivs();
    });
    remakeDivs();
    
        
    //Make the helpbox draggable.
    $("#chatbox").draggable();
}


//Start-up (document.ready) functions run by the user who is trying to solve
//synthesis problems, only.
function clientSetup(typeableReagents) {
    
    //Autocomplete box for reagents. Borrowed from jquery.
    $( "#reagentTyperBox" ).bind( "keydown", function( event ) {
        if ( event.keyCode === $.ui.keyCode.TAB && $( this ).data( "autocomplete" ).menu.active ) {
            event.preventDefault();
        }
    }).autocomplete({
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

    

    
    $("#input").keydown(function(e){
        //If the user presses the enter key (#13) while on the input box,
        //automatically submit line.
        if (e.keyCode == 13) {
            submitLine();
        }
    });
    
    //For chat functionality.    
    $("#submit").click(function(){
        submitLine();
    });
    
    //Close chat.
    $("#close").click(function(){
        $("#chatbox").css("margin-left", "-9999px");
    });
    
    //Submit reagent when enter is pressed
    $("#reagentTyperBox").keydown(function(e){
        if (e.keyCode == 13) {
            updateReagents();
            $("#reagentTyperBox").val("");
            //Also, get rid of the dropdown menu so that they can drag easily.
            $("#addReagent").focus();
        }
    });
    
    //Make big molecules.
    $("#bigMolecule").each(function(){
        $(this).click(function(){
            $(this).css("left", "-9999px");

        });
    });
    
    updateBigMolecule();
    
    
}

//For making a link display the solution
function solutionDisplay() {
    $.ajax({
        type: "GET",
        url: "/orgo/api/displaySolution/",
        data: {},
        success: function(data) {
            drawAllTheThings(data);
            isSolution = true;
        },
    });
}

//Makes the bigMolecule svg show up when you hover over the corresponding regular molecule.
function updateBigMolecule(){
    $(".molecule, #target").each(function(){
        var bigSelector = "#bigMolecule";
        $(this).click(function(){
            if ($(this).html().length < 5){
                //Sorta hackish way of testing whether anything is in the molecule div
                return
            }
            out=$(this).html().replace('height="200px"', 'height="400px"').replace('width="200px"', 'width="400px"');
            image = "<img src=\"https://felixsun.scripts.mit.edu/orgo/static/x.png\" style=\"position:absolute;\"></img>";
            $(bigSelector).css('left', '400px');
            out = image + out;
            $(bigSelector).html(out);
            moleculeVisible = true;

        });
    });
}

var updateReagents = function() {
    //Don't need to send anything back!
    //Update the reagents present in the sidebar with a new reagent.
    //Keep data about that reagent in currentReagents.
    
    //If we're in tutorial mode, advance to the next set of boxes.
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

function submitLine() {
    if (state =! "Chatting") {
        return;
    }
    //Submit a new message to the server.
    $.ajax({
        type: "POST",
        url: "/orgo/chat/helpeechatpoll/",
        data: {'PK': pk,
               'message': $('#input').val()},
        success: update
    });
    //Reset the text input.
    $('#input').val("");
}

function askForHelp() {
    $("#chatbox").css('margin-left', '100px');
    $.ajax({
        type: "GET",
        url: "/orgo/chat/askforhelp/",
        success: function(data) {
            if(data==1){
              $("#helpbox").html("1 person waiting for help, including you");
            }else{
              $("#helpbox").html(data+" people waiting for help, including you");
            }
        }
    });
    setTimeout(keepPolling, REFRESH);
}

function keepPolling() {
    $.ajax({
        type: "GET",
        url: "/orgo/chat/helpeepoll/",
        success: function(data) {
            //$("#allOut").html(data);
            dataObject = jQuery.parseJSON(data);
            if (dataObject.success){
                $("#helpbox").html("Success!  Helped by "+dataObject.helper+"<br / >");
                state = "Chatting";
                pk = dataObject.chatPK;
                getChat();
            } else {
                if(data==1){
                  $("#helpbox").html("1 person waiting for help, including you");
                }else{  
                  $("#helpbox").html(dataObject.queueSize+" people waiting for help, including you");
                }
                setTimeout(keepPolling, REFRESH);
            }

        }
    });
}

function getChat() {
    //Gets chat messages from the server.
    $.ajax({
        type: "POST",
        url: "/orgo/chat/helpeechatpoll/",
        data: {'PK': pk},
        success: update
    });
    setTimeout(getChat, REFRESH);
}

function update(data) {
    //$("#allOut").html(data);
    jsonObject = $.parseJSON(data);
    if (!jsonObject['open']) {
        $("#helpbox").html("Session dropped.  Did the other user log out?");
        state = "Normal";
        setTimeout(function(){
            $("#chatbox").css('margin-left', '-9999px');
            }, 3000);
        return;
    }
    for (i=0; i<jsonObject['length']; i++) {
        $("#helpbox").append(jsonObject[i] + "<br / >");
    }
    //Scroll to bottom.
    if (jsonObject['length'] != 0) {
        $("#helpbox").scrollTop($("#helpbox").prop("scrollHeight")-$("#helpBox").height());
    }
}

function saveProblem() {
    $.ajax({
        type: "GET",
        url: "/orgo/api/saveProblem",
        success: function(data){
            out="Problem saved. <br /> ID: "+data;
            $("#messageArea").html(out);
            $("#messageBox").css("display", "block");
        }
    });
}

function remakeDivs() {
    $("#leftbar").height($(window).height()-60);
    $("#rightbar").height($(window).height()-60);
}

function closeTutorial() {
    tutorial = 0;
    $(".wave3").css("display", "none");
}




