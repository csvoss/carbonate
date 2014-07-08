//Take in a JSON object containing a list of pairs (idnumber, svg representations), and a list of pairs (idnumber1, idnumber2, reagentstext)
//Replaces the innerhtml of anything with id #flowchartDiagram


// jsPlumb.importDefaults({
    // Connector:"StateMachine",
    // PaintStyle:{ lineWidth:3, strokeStyle:"#000000"},
    // ConnectionOverlays : [
        // [ "Arrow", { 
            // location:1,
            // id:"arrow",
            // length:10,
            // foldback:0.5
        // } ],
    // ]
// });
  
// NOTE here we are just using getSelector so we don't have to rewrite the code for each of the supported libraries.
// you can just use the approriate selector from the library you're using, if you want to. like $(".shape) on jquery, for example.
//var shapes = $(".moleculeShape");
    
// loop through them and connect each one to each other one.
// for (var i = 0; i < shapes.length; i++) {
    // for (var j = i + 1; j < shapes.length; j++) {						
        // var conn = jsPlumb.connect({
            // source:shapes[i],  // just pass in the current node in the selector for source 
            // target:shapes[j],
            // parameters:{"reagents":reagents(shapes[i].attr("numberid"), shapes[j].attr("numberid"))}
            // here we supply a different anchor for source and for target, and we get the element's "data-shape"
            // attribute to tell us what shape we should use, as well as, optionally, a rotation value.
            // anchors:[
                // [ "Perimeter", { shape:$(shapes[i]).attr("data-shape"), rotation:$(shapes[i]).attr("data-rotation") }],
                // [ "Perimeter", { shape:$(shapes[j]).attr( "data-shape"), rotation:$(shapes[j]).attr("data-rotation") }]
            // ]
        // });				
        // conn.connection.getOverlay("label").setLabel(conn.getParameters().reagents);
    // }	
// }

