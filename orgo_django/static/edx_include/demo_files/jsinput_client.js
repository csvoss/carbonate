var jsinput_client = (function () {

    // return JSON representation to be used by server-side grader
    function getGrade() {
        return right_answer_for_edx();
    }

    // return JSON representation of persistent state
    function getState() {
        return student_answer_for_edx();
    }

    // process incoming state from jsinput framework
    // This function will be called with 1 argument when JSChannel is not used,
    // 2 otherwise. In the latter case, the first argument is a transaction 
    // object that will not be used here (see http://mozilla.github.io/jschannel/docs/)
    function setState() {
        var stateStr = arguments.length === 1 ? arguments[0] : arguments[1];
        reset_student_answer(stateStr);
    }

    // set up editor inside of div's with class "jade"
    function setup() {
        // Establish a channel only if this application is embedded in an iframe.
        // This will let the parent window communicate with this application using
        // RPC and bypass SOP restrictions.
        if (window.parent !== window) {
            channel = Channel.build({
                window: window.parent,
                origin: "*",
                scope: "JSInput"
            });

            channel.bind("getGrade", getGrade);
            channel.bind("getState", getState);
            channel.bind("setState", setState);
        }
    }

    //////////////////////////////////////////////////////////////////////
    //
    // Module exports
    //
    //////////////////////////////////////////////////////////////////////

    return {
        setup: setup,   // called to initialize jade editors on this page

        // communication to/from edX jsinput framework
        getState: getState,
        setState: setState,
        getGrade: getGrade
    };

}());

// set up communication with parent
$(document).ready(jsinput_client.setup);


