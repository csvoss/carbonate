

var jq = jQuery.noConflict();

jq(document).ready(function() {
    //For making molecules and reactions draggable
    //$( ".molecule" ).draggable({helper: "clone", revert:true, revertDuration: 100});
    $( ".reaction" ).draggable({helper: "clone", revert:true, revertDuration: 100});
    $( ".reagent" ).draggable({helper: "clone", revert:true, revertDuration: 100});
    
    $("#startingMolecule").droppable({
        drop: function(event, ui) {
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
                    if (dataObject.success)
                        $("#successbox").html("<div style=\"background-color:#00FF00\"><h2>SUCCESS!</h2></div>");
                    else
                        $("#successbox").html("<div style=\"background-color:#FFFFFF\"><h2>Not quite.</h2></div>");
                },
            });
        }
    });

    
    //
    
    
    //For typing in autocompleted reagents
    //Update this line by parsing the value of [item for sublist in [REAGENTS[x][1] for x in range(len(REAGENTS)+1) if not x==0] for item in sublist]   in synthProblem.py
    var typeableReagents = ['H2', 'Hydrogen', 'PdC', 'Pd/C', 'Pd|C', 'Pd C', 'palladium', 'EtOH', 'Ethanol', 'Ethyl alcohol', 'C2H5OH', 'HF', 'Hydrogen fluoride', 'Hydrofluoric acid', 'HBr', 'Hydrogen bromide', 'Hydrobromic acid', 'HCl', 'Hydrogen chloride', 'Hydrochloric acid', 'HI', 'Hydrogen iodide', 'Hydroiodic acid', 'CH2Cl2', 'Dichloromethane', 'Fluorine', 'F2', 'Bromine', 'Br2', 'Chlorine', 'Cl2', 'Iodine', 'I2', 'ROOR', 'tBuOOtBu', 'Peroxide', 'Tert-butyl peroxide', 'Di-tert-butyl peroxide', 'mCPBA', 'PhCO3H', 'RCO3H', 'H2SO4', 'Sulfuric acid', 'H2O', 'Water', 'HOH', 'H20', 'HgSO4', 'Hg2+', 'Mercury sulfate', 'BH3', 'Borane', 'THF', 'Tetrahydrofuran', 'NaOH', 'Sodium hydroxide', 'Hydroxide', 'OH-', 'H2O2', 'Hydrogen peroxide', 'oso4', 'osmium tetroxide', 'osmium oxide', 'NMO', 'NMMO', 'N-Methylmorpholine N-oxide', 'Acetone', 'Propanone', '(CH3)2CO', 'Ozone', 'O3', 'Dimethyl sulfide', 'Methylthiomethane', 'Me2S', 'Zn', 'Zinc', 'Lindlar', 'Sodium', 'Na', 'NH3', 'Ammonia', 'Sodium amide', 'Sodamide', 'NaNH2', 'Amide', '1', 'equivalent', 'one', 'heat', 'hv', 'light', 'hnu', 'tert-butoxide', 'KOC(CH3)3'];
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
        var reagentString = $("#reagentTyperBox").val();
        
        $.ajax({
            type: "POST",
            url: "/orgo/api/returnReagentHtml/",
            data: {'reagentString': reagentString},
            success: function(data) {
                if (data != "") 
                $("#reagentsHere").prepend(data);
            },
        });
    }
    

};
    
});




var updateReagents;



