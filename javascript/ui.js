
function cosmoToolsUI(){

    const LIN = 0, LOG = 1;

    var cm        = null, 
        myTable   = null,
        chartURI  = null;

    /** What to do on changing the flatness switch. */
    function onFlat(){
        let flat = document.getElementById("flat").checked;
        document.getElementById("omega-de").disabled = flat;
    }

    /** What to do on changing the neutrino switch. */
    function onNeutrino(){
        let hasnu = document.getElementById("has-neutrino").checked;
        document.getElementById("omega-nu").disabled = !hasnu;
        document.getElementById("num-nu").disabled   = !hasnu;
    }

    /** What to do on changing value of matter density. */
    function onOmegaMatter(){
        if (!document.getElementById("flat").checked){
            return;
        }

        let OmegaMatter = document.getElementById("omega-matter").value;
        document.getElementById("omega-de").value = 1 - OmegaMatter;
    }

    /** Select a specific cosmology model. */
    function selectCosmology(){

        let models = {
                        plank18  : { 
                                        h     : 0.6736, 
                                        Om0   : 0.3153, 
                                        Ob0   : 0.0493, 
                                        Ode0  : 0.6947, 
                                        sigma8: 0.8111, 
                                        ns    : 0.9649, 
                                        Tcmb0 : 2.7255, 
                                        flat  : false, 
                                        has_nu: false
                                    },
                        plank15  : {
                                        h     : 0.6790, 
                                        Om0   : 0.3065, 
                                        Ob0   : 0.0483, 
                                        Ode0  : 0.6935, 
                                        sigma8: 0.8154, 
                                        ns    : 0.9681, 
                                        Tcmb0 : 2.7255, 
                                        flat  : false,
                                        has_nu: false
                                    },
                        wmap08   : {
                                        h     : 0.719, 
                                        Om0   : 0.2581, 
                                        Ob0   : 0.0441, 
                                        Ode0  : 0.7420, 
                                        sigma8: 0.796, 
                                        ns    : 0.963, 
                                        Tcmb0 : 2.7255, 
                                        flat  : false,
                                        has_nu: false
                                    },
                        millanium: {
                                        h     : 0.73, 
                                        Om0   : 0.25, 
                                        Ob0   : 0.045, 
                                        sigma8: 0.9, 
                                        ns    : 1.0, 
                                        Tcmb0 : 2.7255, 
                                        flat  : true,
                                        has_nu: false
                                    },
                        custom   : {
                                        flat  : true,
                                        has_nu: false
                                    },
                    };

        let IDs =  {
                        h     : { id: "hubble",       type: "text" , default: ""    },
                        Om0   : { id: "omega-matter", type: "text" , default: ""    },  
                        Ob0   : { id: "omega-baryon", type: "text" , default: ""    },  
                        sigma8: { id: "sigma8",       type: "text" , default: ""    },
                        ns    : { id: "ps-index",     type: "text" , default: ""    },  
                        Tcmb0 : { id: "cmb-temp",     type: "text" , default: ""    },  
                        Ode0  : { id: "omega-de",     type: "text" , default: ""    },  
                        Onu0  : { id: "omega-nu",     type: "text" , default: 0.0   },  
                        Nnu   : { id: "num-nu",       type: "text" , default: 1.0   }, 
                        flat  : { id: "flat",         type: "check", default: true  },
                        has_nu: { id: "has-neutrino", type: "check", default: false },   
                    };

        let model = models[ document.getElementById("cosmo").value ];
        if (model == undefined) {
            model = models.custom;
        }

        var input = {};

        for (let [key, obj] of Object.entries(IDs)){
            let value = model[ key ];
            if (value == undefined){
                value = obj.default;
            }

            let field = document.getElementById( obj.id );
            
            switch (obj.type) {
                case "text":
                    field.value = value;
                    break;

                case "check":
                    field.checked = value;
            
                default:
                    break;
            }

            input[obj.id] = value;

        }

        onFlat();
        onNeutrino();
        onOmegaMatter();
        
    }

    /** Get the parameters from the request */
    function getParameters(){

        let query = window.location.search.replace('?', '');

        if (query.length == 0){
            var alertBox = document.getElementById('alert-msg');
            alertBox.innerHTML = "<strong>Error!</strong> No cosmology model is specified.";
            alertBox.parentElement.style.display = "block";
            return;
        }

        let paramID = {
                        omegaMatter   : 'Om0',
                        omegaBaryon   : 'Ob0',
                        omegaNeutrino : 'Onu0',
                        numberNeutrino: 'Nnu',
                        omegaDE       : 'Ode0',
                        sigma8        : 'sigma8',
                        psIndex       : 'ns',
                        hubble        : 'h',
                        tempCMB       : 'Tcmb0',
                        w0DE          : 'w0',
                        waDE          : 'wa',
                        flatness      : 'flat',
                    }; 

        var param = { flat: false };

        for (let x of query.split('&')){
            
            let [key, value] = x.split('=');
            key              = paramID[key];
            
            if (key === 'flat'){
                param[key] = ( value === 'flat' ) ? true : false ;
            } else {
                value      = parseFloat(value);
                value      = isNaN(value) ? undefined : value;
                param[key] = value;
            }
        }

        try{
            cm = Cosmo.cosmology(param);
        } catch (err){
            var alertBox = document.getElementById('alert-msg');
            alertBox.innerHTML = "<strong>Error!</strong> " + err.message + ".";
            alertBox.parentElement.style.display = "block";
            return;
        }
    }
    
    /** Create a table of specified values */
    function createTable(value, args){

        /** Return a table of power spectrum values given wavenumber */
        function createPowerTable(args){

            const { ka, kb, z, rows, step } = args;
            
            let lnka = Math.log(ka), lnkb = Math.log(kb);
            let dlnk = (rows > 1) ? (lnkb - lnka) / (rows-1) : 0;

            let table = [];
            for (let i = 0, lnk = lnka; i < rows; i++, lnk += dlnk){
                let k    = Math.exp(lnk);
                let o    = cm.matterPowerSpectrum(k, z, true, true);
                let neff = cm.effectiveIndex(k, z, step);

                table.push([k, o.Pk, o.Dk, o.Tk, neff]);
            }

            return table;
        }

        /** Return a table of correlation values given distance */
        function createCorrfuncTable(args){

            const { ra, rb, z, rows } = args;
            
            let lnra = Math.log(ra), lnrb = Math.log(rb);
            let dlnr = (rows > 1) ? (lnrb - lnra) / (rows-1) : 0;

            let table = [];
            for (let i = 0, lnr = lnra; i < rows; i++, lnr += dlnr){
                let r = Math.exp(lnr);

                table.push([r, cm.matterCorrelation(r, z)]);
            }

            return table;
        }

        /** Returns the mass function table */
        function createMassfuncTable(args){

            const { ma, mb, rows, z, overdensity, lbcorrect } = args;

            let lnma = Math.log(ma), lnmb = Math.log(mb);
            let dlnm = (rows > 1) ? (lnmb - lnma) / (rows-1) : 0.0;

            let table = [];
            for (let i = 0, lnm = lnma; i < rows; i++, lnm += dlnm){

                let m   = Math.exp(lnm);
                let obj = cm.massFunction(m, z, overdensity, 'full');

                let nu  = cm.peakHeight(obj.sigma, z, "s");
                let lb  = cm.linearBias(nu, z, overdensity, lbcorrect);

                table.push([obj.m, obj.r, obj.dndm, obj.dndlnm, obj.f, obj.sigma, obj.dlnsdlnm, lb]);
            }

            return table;
        }

        /** Returns the growth factor table */
        function createGrowthTable(args){

            const { za, zb, rows } = args;

            let dz = (rows > 1) ? (zb - za) / (rows - 1) : 0.0;

            let dplus0 = cm.dplus(0);
            let table  = [];
            for (let i = 0, z = za; i < rows; i++, z += dz){

                table.push([ z, 1.0/(z + 1), cm.dplus(z) / dplus0, cm.growthRate(z) ]);
            }

            return table;
        }

        let tableStruct = {
                            power   : {
                                        title   : "Matter Power Spectrum",
                                        cols    : 5,
                                        colnames: [
                                                    "Wavenumber",
                                                    "Power Spectrum",
                                                    "Dimenssionless Power",
                                                    "Transfer Function",
                                                    "Effective Index"
                                                  ],
                                        func    : createPowerTable, 
                                        activeX : { col: 0, scale: LOG },
                                        activeY : { col: 1, scale: LOG }
                                      },
                            corrfunc: {
                                        title   : "Matter Correlation Function",
                                        cols    : 2,
                                        colnames: [
                                                    "Distance",
                                                    "Correlation"
                                                  ],
                                        func    : createCorrfuncTable,
                                        activeX : { col: 0, scale: LOG },
                                        activeY : { col: 1, scale: LOG }
                                      },
                            massfunc: {
                                        title   : "Halo Mass-function",
                                        cols    : 8,
                                        colnames: [
                                                    "Mass", 
                                                    "Radius", 
                                                    "dndM", 
                                                    "dndlnM", 
                                                    "f(σ)", 
                                                    "σ", 
                                                    "dlnσdlnM",
                                                    "Linear Bias"
                                                  ],
                                        func    : createMassfuncTable,
                                        activeX : { col: 0, scale: LOG },
                                        activeY : { col: 3, scale: LOG }
                                      },
                            growth  : {
                                        title   : "Linear Growth Factor",
                                        cols    : 4,
                                        colnames: [
                                                    "Redshift",
                                                    "Scale Factor",
                                                    "Linear Growth Factor",
                                                    "Linear Growth Rate"
                                                  ],
                                        func    : createGrowthTable,
                                        activeX : { col: 0, scale: LIN },
                                        activeY : { col: 2, scale: LIN }
                                      },
                          };

        
        let plan = tableStruct[value];
        if (plan == undefined){
            return console.log("Error: " + value);
        } 

        let data = plan.func(args);

        myTable = {
                    title   : plan.title,
                    cols    : plan.colnames,
                    rows    : data,
                    shape   : [data.length, plan.cols],
                  };
        
        setChartOptions(    
                            plan.colnames, 
                            plan.activeX.col, 
                            plan.activeY.col, 
                            plan.activeX.scale, 
                            plan.activeY.scale
                        );
        return;
    }

    /** Show or hide advansed options for power */
    function toggleAdvancedOptions(){

        var content = document.getElementById("adv-opt");
        var icon    = document.getElementById("adv-opt-state");
        var display = content.style.display;

        if (display === "none" || display === ""){
            content.style.display = "block";
            icon.innerHTML        = "&bigtriangleup;";
        } else {
            content.style.display = "none";
            icon.innerHTML        = "&bigtriangledown;"
        }
    }

    /** Display an alert message */
    function raiseAlert(msg){

        var alertBox = document.getElementById('calc-alert-msg');
        alertBox.innerHTML = "<strong>Error!</strong> " + msg + ".";
        alertBox.parentElement.style.display = "block";
        return;
    }

    /** Calculate the power spectrum table */
    function calculatePower(){

        document.getElementById("my-viewer").style.display = "none";

        if (cm == null){
            return raiseAlert("No cosmology model is specified");
        }

        let z = parseFloat(document.getElementById("redshift").value);
        if ( isNaN(z) || z <= -1 ){
            return raiseAlert("Invalid or empty value for redshift");
        }

        let ka = parseFloat(document.getElementById("ka").value);
        if ( isNaN(ka) ){
            return raiseAlert("Invalid or empty value for lower k limit");
        }

        let kb = parseFloat(document.getElementById("kb").value);
        if ( isNaN(kb) ){
            return raiseAlert("Invalid or empty value for upper k limit");
        }

        let N = parseInt(document.getElementById("N").value);
        if ( isNaN(N) ){
            return raiseAlert("Invalid or empty value for number of rows");
        }

        var model = document.getElementById("model").value;
        if ( model === "none" ){
            return raiseAlert("A power spectrum model must be selected");
        }

        // advansed options

        let useExactGrowth = document.getElementById("exact-growth").checked;

        var kaNorm = parseFloat(document.getElementById("ka-integ").value);
        if ( isNaN(kaNorm) ){
            return raiseAlert("Invalid or empty value for lower k integration limit");
        }

        var kbNorm = parseFloat(document.getElementById("kb-integ").value);
        if ( isNaN(kbNorm) ){
            return raiseAlert("Invalid or empty value for upper k integration limit");
        }

        var nNorm = parseInt(document.getElementById("num-integ").value);
        if ( isNaN(nNorm) ){
            return raiseAlert("Invalid or empty value for number of points for integration");
        }

        var rstep = parseFloat(document.getElementById("diff-step").value);
        if ( isNaN(rstep) ){
            return raiseAlert("Invalid or empty value for relative step-size");
        }

        // checking the values

        if (ka <= 0 || kb <= 0 || kaNorm <= 0 || kbNorm <= 0 ){
            return raiseAlert("Wavenumber values must be positive");
        }
        else if (kb <= ka || kbNorm <= kaNorm) {
            return raiseAlert("Lower wavenumber must be less than the upper wavenumber");
        }
        else if (N < 1){
            return raiseAlert("Number of rows must be positive and non-zero");
        }
        else if (nNorm < 3 || nNorm % 2 == 0){
            return raiseAlert("Number of points for integration must be an odd number > 3");
        }
        else if (rstep <= 0) {
            return raiseAlert("Relative step-size must be a non-zero positve number");
        }

        cm.setLinearPowerSpectrumModel(model);
        cm.setupPowerQuadArgs({ka: kaNorm, kb: kbNorm, n: nNorm});
        cm.useExactGrowth(useExactGrowth);
        cm.normalize();

        createTable("power", {ka: ka, kb: kb, z: z, rows: N, step : rstep});

        document.getElementById("my-viewer").style.display    = "block";
        document.getElementById("table-viewer").style.display = "none";
        document.getElementById("chart-viewer").style.display = "none";

        return;
    }

    /** Calculate the correlation function table. */
    function calculateCorrfunc() {

        document.getElementById("my-viewer").style.display = "none";

        if (cm == null){

            // console.log("No cosmology model is specified");
            return raiseAlert("No cosmology model is specified");
        }

        let z = parseFloat(document.getElementById("redshift").value);
        if ( isNaN(z) || z <= -1 ){
            return raiseAlert("Invalid or empty value for redshift");
        }

        let ra = parseFloat(document.getElementById("ra").value);
        if ( isNaN(ra) ){
            return raiseAlert("Invalid or empty value for lower r limit");
        }

        let rb = parseFloat(document.getElementById("rb").value);
        if ( isNaN(rb) ){
            return raiseAlert("Invalid or empty value for upper r limit");
        }

        let N = parseInt(document.getElementById("N").value);
        if ( isNaN(N) ){
            return raiseAlert("Invalid or empty value for number of rows");
        }

        var model = document.getElementById("model").value;
        if ( model === "none" ){
            return raiseAlert("A power spectrum model must be selected");
        }

        // advansed options

        let useExactGrowth = document.getElementById("exact-growth").checked;

        var kaNorm = parseFloat(document.getElementById("ka-integ").value);
        if ( isNaN(kaNorm) ){
            return raiseAlert("Invalid or empty value for lower k integration limit");
        }

        var kbNorm = parseFloat(document.getElementById("kb-integ").value);
        if ( isNaN(kbNorm) ){
            return raiseAlert("Invalid or empty value for upper k integration limit");
        }

        var nNorm = parseInt(document.getElementById("num-integ").value);
        if ( isNaN(nNorm) ){
            return raiseAlert("Invalid or empty value for number of points for integration");
        }

        // checking the values

        if (ra <= 0 || rb <= 0){
            return raiseAlert("Distance values must be positive");
        }
        else if (rb <= ra) {
            return raiseAlert("Lower limit must be less than the upper limit");
        }
        else if (kaNorm <= 0 || kbNorm <= 0 ){
            return raiseAlert("Wavenumber values must be positive");
        }
        else if (kbNorm <= kaNorm) {
            return raiseAlert("Lower wavenumber must be less than the upper wavenumber");
        }
        else if (N < 1){
            return raiseAlert("Number of rows must be positive and non-zero");
        }
        else if (nNorm < 3 || nNorm % 2 == 0){
            return raiseAlert("Number of points for integration must be an odd number > 3");
        }

        // console.log(z, N, ka, kb, model, useExactGrowth, kaNorm, kbNorm, nNorm);

        cm.setLinearPowerSpectrumModel(model);
        cm.setupPowerQuadArgs({ka: kaNorm, kb: kbNorm, n: nNorm});
        cm.useExactGrowth(useExactGrowth);
        cm.normalize();

        createTable("corrfunc", {ra: ra, rb: rb, z: z, rows: N});

        document.getElementById("my-viewer").style.display    = "block";
        document.getElementById("table-viewer").style.display = "none";
        document.getElementById("chart-viewer").style.display = "none";

        return;
    }

    /** Calculate the halo mass function table. */
    function calculateMassfunc(){

        document.getElementById("my-viewer").style.display = "none";

        if (cm == null){

            return raiseAlert("No cosmology model is specified");
        }


        let z = parseFloat(document.getElementById("redshift").value);
        if ( isNaN(z) || z <= -1){
            return raiseAlert("Invalid or empty value for redshift");
        }

        let ma = parseFloat(document.getElementById("ma").value);
        if ( isNaN(ma) ){
            return raiseAlert("Invalid or empty value for lower mass limit");
        }

        let mb = parseFloat(document.getElementById("mb").value);
        if ( isNaN(mb) ){
            return raiseAlert("Invalid or empty value for upper mass limit");
        }

        let N = parseInt(document.getElementById("N").value);
        if ( isNaN(N) ){
            return raiseAlert("Invalid or empty value for number of rows");
        }

        let overdensity = document.getElementById("overdens").value;
        overdensity = cm.parseOverdensity(overdensity, z);
        if ( isNaN(overdensity) ){
            return raiseAlert("Invalid or empty value for overdensity");
        }

        var mfmodel = document.getElementById("model").value;
        if ( mfmodel === "none" ){
            return raiseAlert("A mass-function model must be selected");
        }

        // advansed options

        let useExactGrowth = document.getElementById("exact-growth").checked;

        var kaNorm = parseFloat(document.getElementById("ka-integ").value);
        if ( isNaN(kaNorm) ){
            return raiseAlert("Invalid or empty value for lower k integration limit");
        }

        var kbNorm = parseFloat(document.getElementById("kb-integ").value);
        if ( isNaN(kbNorm) ){
            return raiseAlert("Invalid or empty value for upper k integration limit");
        }

        var nNorm = parseInt(document.getElementById("num-integ").value);
        if ( isNaN(nNorm) ){
            return raiseAlert("Invalid or empty value for number of points for integration");
        }

        var psmodel = document.getElementById("ps-model").value;
        var filt    = document.getElementById("filt").value;

        var lbmodel     = document.getElementById("bias-model").value
        var correctBias = document.getElementById("bias-correct").checked;


        // checking the values

        if (ma <= 0 || mb <= 0 ){
            return raiseAlert("Mass values must be positive");
        }
        else if (kaNorm <= 0 || kbNorm <= 0 ){
            return raiseAlert("Wavenumber values must be positive");
        }
        else if (kbNorm <= kaNorm) {
            return raiseAlert("Lower wavenumber must be less than the upper wavenumber");
        }
        else if (N < 1){
            return raiseAlert("Number of rows must be positive and non-zero");
        }
        else if (nNorm < 3 || nNorm % 2 == 0){
            return raiseAlert("Number of points for integration must be an odd number > 3")
        }

        cm.setLinearPowerSpectrumModel(psmodel);
        cm.setFilter(filt);
        cm.setupPowerQuadArgs({ka: kaNorm, kb: kbNorm, n: nNorm});
        cm.useExactGrowth(useExactGrowth);
        cm.normalize();

        cm.setHaloMassfunctionModel(mfmodel);
        cm.setLinearBiasModel(lbmodel);

        createTable("massfunc", {ma: ma, mb: mb, rows: N, z: z, overdensity: overdensity, lbcorrect: correctBias});

        document.getElementById("my-viewer").style.display    = "block";
        document.getElementById("table-viewer").style.display = "none";
        document.getElementById("chart-viewer").style.display = "none";

        return;

    }

    /** Calculate the linear growth factors table. */
    function calculateGrowth() {

        document.getElementById("my-viewer").style.display = "none";

        if (cm == null){

            // console.log("No cosmology model is specified");
            return raiseAlert("No cosmology model is specified");
        }

        let za = parseFloat(document.getElementById("za").value);
        if ( isNaN(za) ){
            return raiseAlert("Invalid or empty value for lower z limit");
        }

        let zb = parseFloat(document.getElementById("zb").value);
        if ( isNaN(zb) ){
            return raiseAlert("Invalid or empty value for upper z limit");
        }

        let N = parseInt(document.getElementById("N").value);
        if ( isNaN(N) ){
            return raiseAlert("Invalid or empty value for number of rows");
        }

        let useExactGrowth = document.getElementById("exact-growth").checked;


        // checking the values

        if (za + 1 <= 0 || zb + 1 <= 0){
            return raiseAlert("Redshift values must be greater than -1");
        }
        else if (zb <= za) {
            return raiseAlert("Lower limit must be less than the upper limit");
        }
        else if (N < 1){
            return raiseAlert("Number of rows must be positive and non-zero");
        }


        cm.useExactGrowth(useExactGrowth);

        createTable("growth", {za: za, zb: zb, rows: N}),

        document.getElementById("my-viewer").style.display    = "block";
        document.getElementById("table-viewer").style.display = "none";
        document.getElementById("chart-viewer").style.display = "none";

        return;
    }

    /** Set the chart options. */
    function setChartOptions(cols, xcol, ycol, xscale, yscale){

        var opts, rad, i;
        
        opts = "";
        for (let i = 0; i < cols.length; i++){
            opts += "<option value=\"" 
                        + i 
                        + "\" " 
                        + ( (i == xcol) ? "selected" : '' )
                        +  ">" 
                        + cols[i] 
                        + "</option>";
        }
        document.getElementById("chart-xval").innerHTML = opts;

        rad = document.getElementsByName("chartXscale");
        for (let i = 0; i < rad.length; i++ ){
            rad[i].checked = (i == xscale);
        }


        opts = "";
        for (let i = 0; i < cols.length; i++){
            opts += "<option value=\"" 
                        + i 
                        + "\" " 
                        + ( (i == ycol) ? "selected" : '' )
                        +  ">" 
                        + cols[i] 
                        + "</option>";
        }
        document.getElementById("chart-yval").innerHTML = opts;

        rad = document.getElementsByName("chartYscale");
        for (let i = 0; i < rad.length; i++ ){
            rad[i].checked = (i == yscale);
        }

    }

    /** Get the chart options. */
    function getChartOptions(){

        let rad, xcol, ycol, xscale, yscale;

        // x-axis

        xcol = document.getElementById("chart-xval").selectedIndex;
        rad = document.getElementsByName("chartXscale");
        for(let i = 0; i< rad.length; i++){
            if (rad[i].checked){
                xscale = i; 
                break;
            }
        }


        // y-axis

        ycol = document.getElementById("chart-yval").selectedIndex;
        rad = document.getElementsByName("chartYscale");
        for(let i = 0; i< rad.length; i++){
            if (rad[i].checked){
                yscale = i; 
                break;
            }
        }

        return { xcol: xcol, ycol: ycol, xscale: xscale, yscale: yscale };
    }

    /** Create a chart with google charts */
    function createChart(){

        /** Create the chart. */
        function makeChart(){

            const { xcol, ycol, xscale, yscale } = getChartOptions();

            var xtitle = myTable.cols[xcol], 
                ytitle = myTable.cols[ycol];
            var title  = xtitle + " vs " + ytitle;

            var data = [[ xtitle, ytitle ]];
            for (let i = 0; i < myTable.shape[0]; i++){
                data.push([     
                                myTable.rows[i][xcol], myTable.rows[i][ycol] 
                          ]);
            }

            data = google.visualization.arrayToDataTable(data);

            var options = {
                                title: title,
                                hAxis: {
                                            title   : xtitle,
                                            logScale: (xscale == LOG) ? true : false,
                                            format  : 'scientific',
                                    },
                                vAxis: {
                                            title   : ytitle,
                                            logScale: (yscale) ? true : false,
                                            format  : 'scientific',
                                    },
                                legend: 'none',
                        };
            
            var chart1 = new google.visualization.LineChart(document.getElementById('my-chart'));
            google.visualization.events.addListener(chart1, 'ready', function () { chartURI = chart1.getImageURI();});
            chart1.draw(data, options);
        }

        try{

            google.charts.load('current',{packages:['corechart']});
            google.charts.setOnLoadCallback(makeChart);
        }
        catch (err){
            console.log(err);
            return -1;
        }
        return 0;
    }

    /** Show the last calculated values as chart */
    function showChart(){

        // var width = (window.innerWidth > 0) ? window.innerWidth : screen.width;
        // var disp  = (width < 800) ? "block" : "flex";

        if (createChart() == -1){
            var alertBox = document.getElementById('alert-msg-viewer');
            alertBox.innerHTML = "<strong>Message!</strong> Failed to create chart!";
            alertBox.parentElement.style.display = "block";
        }

        let btn1 = document.getElementById("chart-btn");
        btn1.className += " active-opt";

        let btn2 = document.getElementById("table-btn");
        btn2.className = btn2.className.replace(" active-opt", "");

        document.getElementById("table-viewer").style.display = "none";
        document.getElementById("chart-viewer").style.display = "block";

        return;
    }

    /** Open the chart as an image on another tab. */
    function downloadImage(){

        if (chartURI == null) {
            
            var alertBox = document.getElementById('alert-msg-viewer');
            alertBox.innerHTML = "<strong>Message!</strong> Failed to create chart!";
            alertBox.parentElement.style.display = "block";
            return;
        }

        // var str = "";
        // str += "<html><head><title>" + myTable.title + "</title>"
        // str += "<style>img{display: block; margin: auto; width: 640px;}</style>"
        // str += "</head><body>";
        // str += "<img src=" + chartURI + " >";
        // str += "</body></html>";

        // return window.open("").document.write(str);

        // get filename
        let fname = document.getElementById("fname-img").value.trim();
        if (fname === ''){
            var alertBox = document.getElementById('alert-msg-viewer');
                alertBox.innerHTML = "<strong>Error!</strong> Filename is not given.";
                alertBox.parentElement.style.display = "block";
                return;
        }
        
        // create a download link
        let a = document.createElement('a');
        a.setAttribute('href', chartURI);
        a.setAttribute('download', fname);
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);

    }

    /** Show the table */
    function showTable(){

        document.getElementById('alert-msg-viewer').parentElement.style.display = "none";

        if (myTable == null){
            return;
        }

        var str = "<tr>";
        for (let colname of myTable.cols){
            str += "<th>" + colname + "</th>";
        }
        str += "</tr>";


        let i = 0;
        for (let row of myTable.rows){

            // if (i === 1000){

            //     var alertBox = document.getElementById('alert-msg-viewer');
            //     alertBox.innerHTML = "<strong>Warning!</strong> Only the first " + i + " lines will be shown. Use the 'View as CSV' option to view full data.";
            //     alertBox.parentElement.style.display = "block";
            //     break;
            // }

            str += "<tr>";
            for (let data of row){
                str += "<td>" + data.toExponential(3) + "</td>";
            }
            str += "</tr>";

            i++;
        }

        document.getElementById("table-1").innerHTML = str;

        let btn1 = document.getElementById("table-btn");
        btn1.className += " active-opt";

        let btn2 = document.getElementById("chart-btn");
        btn2.className = btn2.className.replace(" active-opt", "");

        document.getElementById("table-viewer").style.display = "block";
        document.getElementById("chart-viewer").style.display = "none";
    }

    /** download full data on a new blank page in CSV format */
    function downloadData(){

        if (myTable == null){
            return;
        }

        var str = "";

        str += "#" + myTable.cols.join(", ") + "\n";
        for (let row of myTable.rows){

            let r = [];
            for (let data of row){
                r.push( data.toExponential(8) );
            }
            str += r.join(", ") + "\n";
        }

        // return window.open("").document.write(str);

        // get filename
        let fname = document.getElementById("fname-tbl").value.trim();
        if (fname === ''){
            var alertBox = document.getElementById('alert-msg-viewer');
                alertBox.innerHTML = "<strong>Error!</strong> Filename is not given.";
                alertBox.parentElement.style.display = "block";
                return;
        }
        
        // create a download link
        let a = document.createElement('a');
        a.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(str));
        a.setAttribute('download', fname);
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);

    }

    /** Close viewer tab */
    function closeViewer(){
        document.getElementById("chart-viewer").style.display = "none";
        document.getElementById("table-viewer").style.display = "none";

        let btn1 = document.getElementById("chart-btn");
        btn1.className = btn1.className.replace(" active-opt", "");

        let btn2 = document.getElementById("table-btn");
        btn2.className = btn2.className.replace(" active-opt", "");
    }

    /** Select a calculator mode. */
    function selectMode(cmode){

        let query = window.location.search;
        let url = "./calculator-";

        switch (cmode) {
            case "ps":
                url += "power";
                break;

            case "mf":
                url += "massfunc";
                break;

            case "xi":
                url += "corrfunc";
                break;

            case "zg":
                url += "growth";
                break;
        
            default:
                return;
        }


        url += ".html" + query;
        window.open(url, "_self");
    }
    


    return {
                onFlat               : onFlat,
                onNeutrino           : onNeutrino,
                onOmegaMatter        : onOmegaMatter,
                selectCosmology      : selectCosmology,
                getParameters        : getParameters,
                toggleAdvancedOptions: toggleAdvancedOptions,
                showTable            : showTable,
                showChart            : showChart,
                downloadImage        : downloadImage,
                downloadData         : downloadData,
                closeViewer          : closeViewer,
                calculatePower       : calculatePower,
                calculateCorrfunc    : calculateCorrfunc,
                calculateMassfunc    : calculateMassfunc,
                calculateGrowth      : calculateGrowth,
                selectMode           : selectMode,
            }

}

UI = cosmoToolsUI();