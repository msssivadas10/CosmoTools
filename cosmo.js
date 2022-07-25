/**
 * @file cosmo.js
 * @version 1.1
 * @author 
 * 
 * A javascript library for cosmological calculations.
 * 
 */


/**
 * The cosmology module wrapper.
 */
function cosmologyModule(){

    /** Critical density of the universe in h^2 Msun / Mpc^3 */
    const RHO_CRIT0_ASTRO = 2.77536627E+11;

    /** Critical over-density for spherical collapse */
    const DELTA_C = 1.6864701998411453;

    /** Value of 1/(2 pi^2) */
    const ONE_2SQRPI = 0.05066059182116889;

    /** Speed of light in m/sec */ 
    const C_SI = 2997924580; 

    /** Value of Mpc in meters */
    const MPC = 3.085677581491367e+22;

    /** year (sidereal) in seconds */
    const YEAR = 31558149.8; 

    /** 
     * Error thrown by cosmo module. The error will have a `name` and a `message` attribute.
     */
    function CosmoError(message){

        return {name: 'CosmoError', message: message};
    }


    /**
     * A library of some numerical methods such as integration and interpolation.
     */
    function numerical(){

        function adaptiveSimpsonAux(f, a, b, eps, whole, fa, fb, fm, rec){
            let m = (a + b) * 0.5;
            let h = (b - a) * 0.5;

            let lm = (a + m) * 0.5, rm = (m + b) * 0.5;

            if( (eps/2 == eps) || (a == lm) ){
                return whole;
            }

            let flm = f(lm), frm = f(rm);
            let left  = (h/6) * (fa + 4*flm + fm);
            let right = (h/6) * (fm + 4*frm + fb);
            let delta = left + right - whole

            if (rec <= 0){
                throw CosmoError("Maximum recursion depth reached");
            }

            if (Math.abs(delta) <= 15*eps){
                return left + right + delta / 15.0;
            }
            return adaptiveSimpsonAux(f, a, m, eps/2, left, fa, fm, flm, rec-1) 
                        + adaptiveSimpsonAux(f, m, b, eps/2, right, fm, fb, frm, rec-1); 

        }

        /**
         * Returns the integral of function computed by adaptive Simpson's rule.
         * 
         * @param f Function to integrate.
         * @param a Lower limit of integration.
         * @param b Upper limit of integration.
         * @param eps Accuracy of the result.
         * @param maxrec Maximum recursions to use.
         */
        function adaptiveSimpson(f, a, b, eps, maxrec){
            if (eps == undefined){ eps = 1e-08; }
            if( maxrec == undefined ){ maxrec = 50; }


            let h = b - a;
            if (h == 0){
                return 0.0;
            }

            let fa = f(a), fb = f(b), fm = f(0.5*(a + b));
            let s  = (h/6) * (fa + 4*fm + fb);
            return adaptiveSimpsonAux(f, a, b, eps, s, fa, fb, fm, maxrec);
        }

        /**
         * Returns the integral of function using Simpson's rule, given a pre-computed table of function 
         * evaluations and the stepsize. The function should be evaluated at equidistant points in the 
         * interval of integration.
         * 
         * @param y Function values.
         * @param dx Stepsize. 
         */
        function simpsonSum(y, dx){
            let n = y.length;
            let m = 0;
            if (n%2 == 0){
                n = n - 1;
                m = 1;
            }

            let sum = 0.0;
            for (let j = 1; j < n; j += 2){
                sum += ( y[j-1] + 4*y[j] + y[j+1] );
            }

            return sum * dx / 3;
        }

        /**
         * Create a natural cublic spline object from the given values x and y. The `eval` function 
         * can be then used to evaluate the spline at any intermediate point.
         */
        function cspline(x, y){

            let n = x.length;
            if ( y.length != n ){
                throw CosmoError("x and y must have same length");
            }
            n--;

            let b = Array.apply(null, Array(n)).map(function() {});
            let d = Array.apply(null, Array(n)).map(function() {});
            let c = Array.apply(null, Array(n+1)).map(function() {});


            let h = [];
            for (let i = 0; i < n; i++){
                h.push( x[i+1] - x[i] );
            }

            let alpha = [null];
            for (let i = 1; i < n; i++){
                alpha.push( 3/h[i] * (y[i+1] - y[i]) - 3/h[i-1] * (y[i] - y[i-1]) )
            }

            let l = [0.0], mu = [0.0], z = [0.0];
            for (let i = 1; i < n; i++){
                l.push( 2*(x[i+1] - x[i-1]) - h[i-1] * mu[i-1] ) ;
                mu.push( h[i] / l[i] );
                z.push( (alpha[i] - h[i-1]*z[i-1]) / l[i] ) ;
            }
            l.push(1.0);
            z.push(0.0);
            c[n] = 0.0;

            for (let j = n-1; j > -1; j--){
                c[j] = z[j] - mu[j] * c[j+1];
                b[j] = ( y[j+1] - y[j] )/h[j] - (c[j+1] + 2*c[j]) * h[j]/3;
                d[j] = (c[j+1] - c[j]) / (3*h[j]);
            }

            var s = [];
            for (let i = 0; i < n; i++){
                s.push( {x: x[i], y: y[i], b: b[i], c: c[i], d: d[i]} );
            }


            var bounds = [ x[0], x[n] ];


            function eval(x){

                if (x < bounds[0] || x > bounds[1]){
                    throw CosmoError("x value out of bounds");
                }

                for (let i = 0; i < s.length; i++){
                    let xa = s[i].x;
                    
                    let xb = bounds[1];
                    if (i < s.length-1){ 
                        xb = s[i+1].x; 
                    }

                    if (xa <= x && x <= xb){
                        let a = s[i].y, b = s[i].b, c = s[i].c, d = s[i].d;
                        let t = x - xa;

                        return a + b*t + c*t*t + d*t*t*t;
                    }
                }
            }

            var publicAPI = {
                                eval: eval,
                            }

            return publicAPI;
        }


        var publicAPI = {
                            simpson: adaptiveSimpson,
                            simpsonSum: simpsonSum,
                            cspline: cspline,
                        };

        return publicAPI;
    }

    /**
     * A library of some available linear matter power spectrum/transfer function models. Available 
     * models are Eisenstein & Hu (1998) models with and without baryon oscillations, Bardeen et al
     * (BBKS) with Sugiyama correction (1995). To set an active power spectrum model use `setActive` 
     * function with a model key (`eisenstein98_bao`, `eisenstein98_zb`, `sugiyama95` respectively).
     */
    function power_spectrum(cm){

        var parameters = cm.parameters;

        var Om0  = parameters.Om0, 
            Ob0  = parameters.Ob0, 
            Onu0 = parameters.Onu0, 
            Nnu  = parameters.Nnu, 
            h    = parameters.h, 
            theta = parameters.Tcmb0 / 2.7;

        var h2   = h*h,
            Omh2 = Om0*h2,
            Obh2 = Ob0*h2;
        var wt_b = Ob0 / Om0;  // fraction of baryons
        var wt_c  = 1 - wt_b;  // fraction of cold dark matter


        // parameters for power spectrum models
        
        var active; // active / selected power spectrum model

        var gamma_eff;      // for sugiyam95
        var s, alpha_gamma; // for eisenstein98
        var k_eq, k_silk, alpha_b, alpha_c, beta_b, beta_c, beta_node; // for eisenstein98_bao


        function sugiyama95_computeParam(){

            gamma_eff = (Om0*h) / Math.exp( Ob0 + Math.sqrt(2*h) * Ob0 / Om0 );
        }

        function sugiyama95(k, z){

            if (k < 1E-09){
                return 1.0;
            }

            let q = k / gamma_eff;

            return Math.log(1 + 2.34*q) / (2.34*q) * Math.pow( 1 + 3.89*q + Math.pow(16.1*q, 2) + Math.pow(5.46*q, 3) + Math.pow(6.71*q, 4), -0.25 );
        }

        function eisenstein98_zb_computeParam(){

            // sound horizon : eqn. 26
            s = 44.5*Math.log(9.83 / Omh2) / Math.sqrt(1 + 10*Math.pow(Obh2, 3/4));

            // eqn. 31
            alpha_gamma = 1 - 0.328*Math.log(431*Omh2)*wt_b + 0.38*Math.log(22.3*Omh2)*wt_b*wt_b;
        }

        function eisenstein98_zb(k, z){

            // eqn. 30
            let gamma_eff = Om0*h*(alpha_gamma + (1 - alpha_gamma)/(1 + Math.pow(0.43*k*s, 4)));

            // eqn. 28
            let q = k*(theta*theta/ gamma_eff);

            // eqn. 29
            let l_0 = Math.log(2*Math.E + 1.8*q);
            let c_0 = 14.2 + 731.0 / (1 + 62.5*q);
            return l_0 / (l_0 + c_0*q*q);
        }

        function eisenstein98_bao_computeParam(){

            // redshift at equality : eqn. 2 (corrected)
            let zp1_eq = (2.50e+04)*Omh2 / Math.pow(theta, 4);

            // wavenumber at equality : eqn. 3
            k_eq = (7.46e-02)*Omh2 / Math.pow(theta, 2);

            // redshift at drag epoch : eqn 4
            let c1  = 0.313*(1 + 0.607*Math.pow(Omh2,0.674)) / Math.pow(Omh2, 0.419);
            let c2  = 0.238*Math.pow(Omh2, 0.223);
            let z_d = 1291.0*Math.pow(Omh2, 0.251)*(1 + c1*Math.pow(Obh2, c2)) / (1 + 0.659*Math.pow(Omh2, 0.828));

            // baryon - photon momentum density ratio : eqn. 5
            let R_const = 31.5*(Obh2 / Math.pow(theta, 4)) * 1000;
            let R_eq    = R_const / zp1_eq;     // at equality epoch
            let R_d     = R_const / (1 + z_d);  // at drag epoch

            // sound horizon : eqn. 6
            s = (2/3/k_eq)*Math.sqrt(6/R_eq)*Math.log((Math.sqrt(1 + R_d) + Math.sqrt(R_eq + R_d)) / (1 + Math.sqrt(R_eq)));

            // silk scale : eqn. 7
            k_silk = 1.6*Math.pow(Obh2, 0.52)*Math.pow(Omh2, 0.73)*(1 + Math.pow(10.4*Omh2, -0.95));

            // eqn. 11
            let a1  = (1 + Math.pow(32.1*Omh2, -0.532))*Math.pow(46.9*Omh2, 0.670);
            let a2  = (1 + Math.pow(45.0*Omh2, -0.582))*Math.pow(12.0*Omh2, 0.424);
            alpha_c = Math.pow(a1, -wt_b) * Math.pow(a2, -Math.pow(wt_b, 3));

            // eqn. 12
            let b1  = 0.944 / (1 + Math.pow(458.0*Omh2, -0.708));
            let b2  = Math.pow(0.395*Omh2, -0.0266);
            beta_c  = 1 / (1 + b1*(Math.pow(wt_c, b2) - 1));

            // eqn. 15
            let y   = zp1_eq / (1 + z_d);
            let y1  = Math.sqrt(1 + y);
            let Gy  = y*( -6*y1 + (2 + 3*y) * Math.log((y1 + 1) / (y1 - 1)) );

            // eqn. 14
            alpha_b = 2.07*(k_eq*s)*Gy*Math.pow(1 + R_d, -3.0/4);

            // eqn. 24
            beta_b  = 0.5 + wt_b + (3 - 2*wt_b)*Math.sqrt(Math.pow(17.2*Omh2, 2) + 1);

            // eqn. 23
            beta_node = 8.41*Math.pow(Omh2, 0.435);
        }

        function eisenstein98_bao(k, z){
            
            k = k*h; // convert wavenumber from h/Mpc to 1/Mpc
            
            let q = k/(13.41*k_eq)  // eqn. 10
            let x = k*s             // new variable

            // eqn. 18
            let f = 1 / (1 + (x/5.4)**4);

            // eqn. 19 and 20
            let l_beta = Math.log(Math.E + 1.8*beta_c*q);

            let c_no_alpha = 14.2 + 386.0 / (1 + 69.9*Math.pow(q, 1.08));
            let t_no_alpha = l_beta / (l_beta + c_no_alpha*q*q);

            let c_alpha    = 14.2 / alpha_c + 386.0 / (1 + 69.9*Math.pow(q, 1.08));
            let t_alpha    = l_beta / (l_beta + c_alpha*q*q);

            // cold-dark matter part : eqn. 17
            let tc = f*t_no_alpha + (1 - f)*t_alpha;

            // eqn. 22
            let s_tilde   = s / Math.pow(1 + (beta_node / x)**3, 1.0/3);
            let x_tilde   = k*s_tilde;

            // eqn. 19 and 20 again
            let l_no_beta = Math.log(Math.E + 1.8*q);
            let t_nothing = l_no_beta / (l_no_beta + c_no_alpha*q*q);

            // baryonic part : eqn. 21
            let j0 = Math.sin(x_tilde) / x_tilde; // zero order spherical bessel
            let tb = (t_nothing / (1 + Math.pow(x / 5.2, 2)) + ( alpha_b / (1 + Math.pow(beta_b / x, 3))) * Math.exp(-Math.pow(k / k_silk, 1.4))) * j0;

            return wt_b * tb + wt_c * tc; // full transfer function : eqn. 16
        }

        /**
         * Set the active power spectrum model. Value for the model argument muat be any of 
         * `eisenstein98_bao`, `eisenstein98_zb` or `sugiyama95`.
         */
        function setActive(model){

            switch (model) {
                case 'sugiyama95':
                    sugiyama95_computeParam();
                    break;

                case 'eisenstein98_zb':
                    eisenstein98_zb_computeParam();
                    break;

                case 'eisenstein98_bao':
                    eisenstein98_bao_computeParam();
                    break;
            
                default:
                    throw CosmoError("unknown power spectrum model: " + model);
                    break;
            }

            active = model;
        }

        /**
         * Returns the transfer function using the active model.
         * 
         * @param k Wavenumber in h/Mpc.
         * @param z Redshift. 
         */
        function transfer(k, z){

            if (active == undefined){
                throw CosmoError("active power spectrum model is not set");
            }

            switch (active) {
                case 'sugiyama95':
                    return sugiyama95(k, z);

                case 'eisenstein98_zb':
                    return eisenstein98_zb(k, z);

                case 'eisenstein98_bao':
                    return eisenstein98_bao(k, z);
            
                default:
                    throw CosmoError("unknown power spectrum model: " + model);
                    break;
            }

        }

        var publicAPI = {
                            setActive: setActive,
                            transfer : transfer,
                        }

        return publicAPI;

    }

    /**
     * Returns the value of a filter function to smooth the density field. The filter function can be 
     * selected by the `model` parameter, whose value should be any of `tophat` or `gauss` for a 
     * spherical tophat or Gaussian filter.
     * 
     * @param x Input argument to the filter function.
     * @param j Order of the derivative. 0 means no deriavative.
     * @param model Filter function to use.
     */
    function filter(x, j, model){

        /** Spherical tophat filter function */
        function tophat(x, j){
            if (j == undefined){ j = 0; }

            switch (j) {
                case 1:
                    return 3*( (x*x - 3)*Math.sin(x) + 3*x*Math.cos(x) ) / (x*x*x*x);
            
                default:
                    return 3*( Math.sin(x) - x*Math.cos(x) ) / (x*x*x);
            }
        }

        /** gaussian filter function */
        function gauss(x, j){
            if (j == undefined){ j = false; }

            switch (j) {
                case 1:
                    return Math.exp(0.5*x*x) * x;
            
                default:
                    return Math.exp(0.5*x*x);
            }
        }    

        if (model == undefined){
            model = 'tophat';
        }

        switch (model) {
            case 'gauss':
                return gauss(x, j);
        
            default:
                return tophat(x, j);;
        }
    }

    /**
     * Returns a mass-function generator object given a cosmology model. The returned object can be 
     * set to compute mass-function based on one of the pre-defined models. The `massFunction` method 
     * can be then used to compute the mass-function values.
     */
    function mass_function(cm){

        var cosmo = cm; // cosmology model

        var num = numerical(); // for cubic splines

        var tinker08_A  = num.cspline(
                                [200,   300,   400,   600,   800,   1200,  1600,  2400,  3200 ],
                                [0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260],
                            );
        var tinker08_a  = num.cspline(
                                [200,   300,   400,   600,   800,   1200,  1600,  2400,  3200 ],
                                [1.47,  1.52,  1.56,  1.61,  1.87,  2.13,  2.30,  2.53,  2.66 ],
                            );
        var tinker08_b  = num.cspline(
                                [200,   300,   400,   600,   800,   1200,  1600,  2400,  3200 ],
                                [2.57,  2.25,  2.05,  1.87,  1.59,  1.51,  1.46,  1.44,  1.41 ],
                            );
        var tinker08_c  = num.cspline(
                                [200,   300,   400,   600,   800,   1200,  1600,  2400,  3200 ],
                                [1.19,  1.27,  1.34,  1.45,  1.58,  1.80,  1.97,  2.24,  2.44 ],
                            );

        var model = 'tinker08';

        const ALL_MODELS = ['press74', 'sheth01', 'jenkins01', 'reed03', 'warren06', 'reed07', 'tinker08', 'crocce10', 'courtin10'];

        /** Mass-function by Press & Schechter (1974). */
        function press74(sigma, z, delta){

            let nu = DELTA_C / sigma;
            let f  = Math.sqrt( 2 / Math.PI ) * nu * Math.exp( -0.5 * nu * nu );
            return f;
        }

        /** Mass-function by Sheth et al (2001). */
        function sheth01(sigma, z, delta){

            let A = 0.3222, a = 0.707, p = 0.3;

            let nu = DELTA_C / sigma;
            let f  = A * Math.sqrt( 2*a / Math.PI ) * nu * Math.exp( -0.5 * a * nu * nu ) * ( 1.0 + Math.pow( nu*nu / a, -p ) );
            return f;
        }

        /** Mass-function Jenkins et al (2001). */
        function jenkins01(sigma, z, delta){

            let f = 0.315 * Math.pow( -Math.abs( Math.log( 1/sigma ) + 0.61 ), 3.8 );
            return f;
        }

        /** Mass-function by Reed et al (2003). */
        function reed03(sigma, z, delta){

            let fst = sheth01(sigma, z, delta);
            let f   = fst * Math.exp( -0.7 / ( sigma * Math.pow( Math.cosh( 2*sigma ), 5) ) );
            return f;
        }

        /** Mass-function by Warren et al (2006). */
        function warren06(sigma, z, delta){

            let A = 0.7234, a = 1.625, b = 0.2538, c = 1.1982;

            let f = A * ( Math.pow(sigma, -a) + b ) * Math.exp( -c / (sigma * sigma) );
            return f;
        }

        /** Mass-function by Reed et al (2007). */
        function reed07(sigma, z, delta){

            if (args.z == undefined){
                throw CosmoError("reed07 mass function is redshift dependent");
            }

            let A = 0.310, c = 1.08, ca = 0.764, p = 0.3;

            let omega = Math.sqrt( ca ) * DELTA_C / sigma;
            let G1    = Math.exp( -0.5*( Math.log( omega ) - 0.788 )**2 / 0.6**2 );
            let G2    = Math.exp( -0.5*( Math.log( omega ) - 1.138 )**2 / 0.2**2 );

            let r     = cosmo.radius( sigma, z );
            let neff  = -2.0 * cosmo.dlnsdlnr( r, 0.0, args ) - 3.0;
            
            // eqn. 12
            let f = A * omega * Math.sqrt( 2.0/Math.PI ) 
                        * Math.exp( -0.5*omega - 0.0325*omega**p / ( neff + 3 )**2 )
                        * ( 1.0 + 1.02*omega**( 2*p ) + 0.6*G1 + 0.4*G2 );
            return f;
        }

        /** Mass-function by Tinker et al (2008). */
        function tinker08(sigma, z, delta){

            if (z == undefined){
                z = 0.0;
            }
            let zp1 = z + 1;

            if (delta == undefined){
                delta = 200;
            }

            if (delta < 200 || delta > 3200){
                throw CosmoError("overdensity must be within 200 and 3200");
            }

            let A = tinker08_A.eval(delta), 
                a = tinker08_a.eval(delta), 
                b = tinker08_b.eval(delta), 
                c = tinker08_c.eval(delta);

            let alpha = Math.pow( 10.0, -Math.pow( 0.75 / Math.log10( delta / 75 ), 1.2 ) ); // eqn 8    

            A = A * Math.pow(zp1, -0.14);  // eqn 5
            a = a * Math.pow(zp1, -0.06);  // eqn 6  
            b = b * Math.pow(zp1, -alpha); // eqn 7 

            let f = A * ( 1 + Math.pow( b / sigma, a ) ) * Math.exp( -c / sigma**2 ); // eqn 3
            return f;
        }

        /** Mass-function by Crocce et al (2010). */
        function crocce10(sigma, z, delta){

            if (z == undefined){
                z = 0.0;
            }
            let zp1 = args.z + 1;

            let Az = 0.580 * Math.pow( zp1, -0.130 );
            let az = 1.370 * Math.pow( zp1, -0.150 );
            let bz = 0.300 * Math.pow( zp1, -0.084 );
            let cz = 1.036 * Math.pow( zp1, -0.024 );
            
            let f = Az * ( Math.pow(sigma, -az) + bz ) * Math.exp( -cz / (sigma*sigma) );
            return f;
        }

        /** Mass-function by Courtin et al (2010). */
        function courtin10(sigma, z, delta){

            let A = 0.348, a = 0.695, p = 0.1;

            let nu = DELTA_C / sigma;
            let f  = A * Math.sqrt( 2*a / Math.PI ) * nu * Math.exp( -0.5 * a * nu * nu ) * ( 1.0 + Math.pow( nu*nu / a, -p ) );
            return f;
        }

        var models = {
                            press74  : press74  ,
                            sheth01  : sheth01  ,
                            jenkins01: jenkins01,
                            reed03   : reed03   ,
                            warren06 : warren06 ,
                            reed07   : reed07   ,
                            tinker08 : tinker08 ,
                            crocce10 : crocce10 ,
                            courtin10: courtin10,
                    }

        /**
         * Set the active halo mass-function model to use. 
         */
        function setModel(m){

            let available = false;
            for (let mi of ALL_MODELS){
                if (m == mi){
                    
                    available = true;
                    break;
                }
            }

            if (!available){
                throw CosmoError("invalid mass function model name: " + m);
                return;
            }

            model = m;
        }

        /**
         * Returns the fitting function based on the active model set.
         * 
         * @param sigma Variance value.
         * @param z Redshift
         * @param delta Overdensity value with respect to the mean background density. 
         */
        function fit(sigma, z, delta){

            // var modelObj = models[model];

            // if ( modelObj == undefined ){
            //     throw CosmoError("invalid mass function model name: " + model);
            // }

            // if ( (modelObj.z_depend == true) && (z == undefined) ){
            //     throw CosmoError("mass function requires redshift value");
            // }
            // if ( (modelObj.mdef == 'so') && (delta == undefined) ){
            //     throw CosmoError("mass function requires an overdensity value");
            // } 

            // return modelObj.f(sigma, z, delta);

            switch (model) {
                case 'press74':
                    return press74(sigma, z, delta);
                    
                case 'sheth01':
                    return sheth01(sigma, z, delta);

                case 'jenkins01':
                    return jenkins01(sigma, z, delta);

                case 'reed03':
                    return reed03(sigma, z, delta);

                case 'warren06':
                    return warren06(sigma, z, delta);

                case 'reed07':
                    return reed07(sigma, z, delta);

                case 'tinker08':
                    return tinker08(sigma, z, delta);

                case 'crocce10':
                    return crocce10(sigma, z, delta);

                case 'courtin10':
                    return courtin10(sigma, z, delta);
            
                default:
                    throw CosmoError("invalid mass function model name: " + model);
            }
        }

        /**
         * Returns the halo mass-function in the format specified by the `fmt` argument, or an object 
         * containing the fields `m` (mass in Msun/h), `r` (Lagrangian radius in Mpc/h), `dndm`, `dndlnm` 
         * `f` (fitting function), `s` (variance/sigma) and `dlnsdlnm` (logarithmic derivative of the 
         * variance w.r.to mass).
         * 
         * @param m Mass of the halo in Msun/h.
         * @param z Redshift.
         * @param overdensity Overdensity value w.r.to mean background density. 
         * @param fmt Output format - any of `dndm`, `dndlnm`, `f` or `full`.
         */
        function massFunction(m, z, overdensity, fmt){

            if (z == undefined){ z = 0.0; }

            if (overdensity == undefined){ overdensity = 200; }

            if (fmt == undefined){ fmt = 'dndlnm'; }


            let r = cosmo.lagrangianR(m);
            let s = Math.sqrt( cosmo.variance(r, z) );
            let dlnsdlnm = cosmo.dlnsdlnr(r, z) / 3;

            let f      = fit(s, z, overdensity);
            let dndlnm = f * (cosmo.rho_m(z) / m) * (-dlnsdlnm);
            let dndm   = dndlnm / m;

            switch (fmt) {
                case 'f':
                    return f;

                case 'dndm':
                    return dndm;

                case 'dndlnm':
                    return dndlnm;

                case 'full':
                    return {
                                m: m, r: r, sigma: s, dlnsdlnm: dlnsdlnm,
                                f: f, dndm: dndm, dndlnm: dndlnm,
                            };
            
                default:
                    throw CosmoError("invalid value for fmt: " + fmt);;
            }
        }

        var publicAPI = {
                            models: models,
                            setModel: setModel,
                            fit: fit,
                            massFunction: massFunction,
                        };

        return publicAPI;
    }

    /**
     * Returns a linear bias function generator object given a cosmology model. The returned object can be 
     * set to compute halo bias values based on one of the pre-defined models. 
     */
    function linear_bias(cm){

        var cosmo = cm;

        /**
         * Linear bias model by Cole & Kaiser (1989), Mo & White (1996).
         */
        function cole89(nu, z, delta, correct){

            return 1. + (nu**2 - 1.) / DELTA_C;
        }

        /**
         * Linear bias model by Sheth et al (2001).
         */
        function sheth01(nu, z, delta, correct){

            let a  = 0.707,
                b  = 0.5,
                c  = 0.6;
            let sqrt_a = Math.sqrt(a);
            let anu2   = a * nu**2;

            return 1.0 + 1.0 / sqrt_a / DELTA_C
                    * ( 
                        sqrt_a * anu2 
                        + sqrt_a * b * anu2**( 1-c ) 
                        - anu2**c / ( anu2**c + b * ( 1-c ) * ( 1-0.5*c ) )
                      )
        }

        /**
         * Linear bias model by Jing (1998).
         */
        function jing98(nu, z, delta, correct){

            let r    = cosmo.radius(DELTA_C / nu, z);
            let neff = cosmo.effectiveIndex(2 * Math.PI / r, z)

            return Math.pow( 0.5 / nu**4 + 1, 0.06 - 0.02*neff ) * ( 1 + ( nu**2 - 1 ) / DELTA_C );
        }

        /**
         * Linear bias model by Seljak et al (2004).
         */
        function seljak04(nu, z, delta, correct){

            if (correct == undefined){ correct = false; }

            let r     = cosmo.radius(nu / DELTA_C, z);
            let rstar = cosmo.radius(1., z);
            let x     = (r / rstar)**3;

            let b = 0.53 + 0.39*x**0.45 + 0.13 / ( 40.0*x + 1 ) + 5E-4*x**1.5;

            if (correct){

                let Om0    = cosmo.parameters.Om0,
                    h      = cosmo.parameters.h,
                    ns     = cosmo.parameters.ns,
                    sigma8 = cosmo.parameters.sigma8;
                
                let alpha_s = 0.0; 

                b += Math.log10(x) * (
                                        0.4*( Om0 - 0.3 + ns - 1.0 ) 
                                         + 0.3*( sigma8 - 0.9 + h - 0.7 ) 
                                         + 0.8*alpha_s
                                    );
            }

            return b;
        }

        /**
         * Linear bias model given by Tinker et al. (2010).
         */
        function tinker10(nu, z, delta, correct){

            let y = Math.log10(delta);
            let x = Math.exp( -Math.pow(4./y, 4) );
            let A  = 1.0 + 0.24 * y * x,
                a  = 0.44 * y - 0.88,
                B  = 0.183,
                b  = 1.5,
                C  = 0.019 + 0.107 * y + 0.19 * x,
                c  = 2.4;

            return 1.0 - A * nu**a / ( nu**a + DELTA_C**a ) + B * nu**b + C * nu**c;
        }

        const models = { 
                            cole89  : cole89, 
                            sheth01 : sheth01, 
                            jing98  : jing98,
                            seljak04: seljak04,
                            tinker10: tinker10,
                       };

        var model = 'tinker10';

        const ALL_MODELS = [ 'cole89', 'sheth01', 'jing98', 'seljak04', 'tinker10' ];

        /**
         * Set the active bias model to use. 
         */
        function setActive(m){

            let available = false;
            for (let mi of ALL_MODELS){

                if (mi == m){
                    available = true;
                    break;
                }
            }

            if (!available){
                throw CosmoError("invalid bias model: " + m);
            }

            model = m;
            return;
        }

        /**
         * Returns the linear bias using the active model set.
         * @param nu Input argument, equal to $\sigma / \delta_c$
         * @param z  Redshift
         * @param overdensity Overdensity value w.r.to the mean background density.
         * @param correct Whether to use correction. Only for model `seljak04`. 
         */
        function bias(nu, z, overdensity, correct){

            if (z == undefined){ 
                z = 0.0; 
            }

            if (overdensity == undefined){ 
                overdensity = 200; 
            }

            if (correct == undefined){
                correct = false;
            }

            switch (model) {
                case 'cole89':
                    return cole89(nu, z, overdensity, correct);

                case 'sheth01':
                    return sheth01(nu, z, overdensity, correct);

                case 'jing98':
                    return jing98(nu, z, overdensity, correct);

                case 'seljak04':
                    return seljak04(nu, z, overdensity, correct);

                case 'tinker10':
                    return tinker10(nu, z, overdensity, correct);
                    
                default:
                    throw CosmoError("invalid bias model: " + model);
            }

        }

        var publicAPI = { 
                            setModel: setActive,
                            bias    : bias,
                            models  : models,
                        }

        return publicAPI;
    }

    /**
     * The cosmology object constructor. This create an object representing a cosmology model and this 
     * object can be used for generating power spectrum and other calculations.
     */
    function cosmology(cm){

        var Om0, Ob0, Ode0, Onu0, Ok0;
        var w0, wa;
        var Mnu, Nnu;
        var Tcmb0;
        var sigma8;
        var ns;
        var flat;
        var h;
        var A = 1.0;


        function init(args){


            if ( args.h == undefined ){
                throw CosmoError("Hubble parameter is a required value");
                return;
            }
        
            h = args.h;
        
            if( h < 0.0 ){
                throw CosmoError("Hubble parameter must be postive");
                return;
            }
        
        
            if ( args.Om0 == undefined ){
                throw CosmoError("Matter density parameter is a required value");
                return;
            }
        
            Om0 = args.Om0;
        
            if( Om0 < 0.0 ){
                throw CosmoError("Matter density parameter must be postive");
                return;
            }
        
        
            if ( args.Ob0 == undefined ){
                throw CosmoError("Baryon density parameter is a required value");
                return;
            }
        
            Ob0 = args.Ob0;
        
            if( Ob0 < 0.0 ){
                throw CosmoError("Baryon density parameter must be postive");
                return;
            }
            else if ( Ob0 > Om0 ){
                throw CosmoError("Baryon density parameter cannot greater than matter density");
                return;
            }
        
        
            Onu0 = 0.0;
            Nnu  = 0.0;
            Mnu  = 0.0;
        
            if ( args.Onu0 != undefined ){
        
                if( Onu0 < 0.0 ){
                    throw CosmoError("Massive neutrino density parameter must be postive");
                    return;
                }
                else if ( Ob0 + Onu0 > Om0 ){
                    throw CosmoError("Sum of baryon and neutrino density parameter cannot greater than matter density");
                    return;
                }
        
                Onu0 = args.Onu0;
        
                if ( args.Nnu == undefined ){
                    throw CosmoError("Number of neutrino species is a required value");
                    return;
                }
        
                Nnu = args.Nnu;
        
                if ( Nnu <= 0.0 ){
                    throw CosmoError("Number of massive neutrinos must be positive");
                    return;
                }
        
                Mnu = 91.5 * Onu0 / Nnu * Math.pow( h, 2 );
        
            }
        
        
            flat = true;
        
            if ( args.flat != undefined ){
                flat = args.flat;
            }
        
            Ode0 = 1 - Om0;
            Ok0  = 0.0;
        
            if ( !flat ){
        
                if ( args.Ode0 == undefined ){
                    throw CosmoError("Dark-energy density parameter is a required value");
                    return 0;
                }
        
                Ode0 = args.Ode0;
                Ok0  = 1 - Om0 - Ode0;
        
            }
        
            if( Ode0 < 0.0 ){
                throw CosmoError("Dark energy density parameter must be postive");
                return;
            }
        
        
            w0 = -1.0;
            wa =  0.0;
        
            if ( args.w0 != undefined ){
                w0 = args.w0;
            }
            if ( args.wa != undefined ){
                wa = args.wa;
            }
        
        
            Tcmb0 = 2.725;
        
            if ( args.Tcmb0 != undefined ){
                Tcmb0 = args.Tcmb0;
            }
            if ( Tcmb0 <= 0.0 ){
                throw CosmoError("CMB temperature must be positive");
                return;
            }
        
        
            if ( args.sigma8 == undefined ){
                throw CosmoError("sigma-8 parameter is a required value");
                return;
            }
            sigma8 = args.sigma8;
        
            if( sigma8 < 0.0 ){
                throw CosmoError("Sigma-8 parameter must be postive");
                return;
            }
        
        
            if ( args.ns == undefined ){
                throw CosmoError("Power spectrum index is a required value");
                return;
            }
            ns = args.ns;
        }

        init(cm);

        var model = { 
                        Om0: Om0, 
                        Ob0: Ob0, 
                        Onu0: Onu0, 
                        Ode0: Ode0, 
                        Ok0: Ok0, 
                        w0: w0, 
                        wa: wa,
                        Nnu: Nnu, 
                        Mnu: Mnu, 
                        ns: ns, 
                        sigma8: sigma8, 
                        h: h, 
                        Tcmb0: Tcmb0, 
                        flat: flat,
                    };

        
        /** Library for muerical methods */
        var numeric = numerical();

        /** Library for power spectrum models */
        var ps = power_spectrum({parameters: model, dplus: dplus});
        ps.setActive('eisenstein98_zb');

        var exactgf = false; // use exact growth factors

        var _filter = 'tophat'; // filter to use

        var kquadargs = {ka: 1e-06, kb: 1e+06, n: 10001};

        var Pktab;
        createPkTable(kquadargs);


        /******************************************/
        /* public API functions                   */
        /******************************************/

        // density calculations

        /**
         * Returns the matter density at redshift z in units of critical density.
         */
        function Om(z){

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;

            let x = Om0 * zp1 * zp1 * zp1;
            let y = x + Ode0 * Math.pow(zp1, 3 + 3*w) + Ok0 * zp1 * zp1;
            return x / y;

        }

        /**
         * Returns the baryon density at redshift z in units if critical density.
         */
        function Ob(z){

            return Om(z) * Ob0 / Om0;
        }

        /**
         * Returns the massive neutrino density at redshift z in units of critical density.
         */
        function Onu(z){

            return Om(z) * Onu0 / Om0;
        }

        /**
         * Returns the dark-energy density at redshift z in units of critical density.
         */
        function Ode(z){

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;

            let x = Ode0 * Math.pow(zp1, 3 + 3*w);
            let y = Om0 * zp1 * zp1 * zp1 + x + Ok0 * zp1 * zp1;
            return x / y;

        }

        /**
         * Returns the curvature density at redshift in units of critical density.
         */
        function Ok(z){

            if (flat){ 
                return 0.0; 
            }

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;

            let x = Ok0 * zp1 * zp1;
            let y = Om0 * zp1 * zp1 * zp1 + Ode0 * Math.pow(zp1, 3 + 3*w) + x;
            return x / y;

        }

        /**
         * Returns the matter density at redshift z in h^2 Msun/Mpc^3. 
         */
        function rho_m(z){

            let zp1 = z + 1;
            return RHO_CRIT0_ASTRO * Om0 * zp1 * zp1 * zp1;

        }

        /**
         * Returns the dark energy density at redshift z in h^2 Msun/Mpc^3.
         */
        function rho_de(z){

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;

            return RHO_CRIT0_ASTRO * Ode0 * Math.pow(zp1, 3 + 3*w)
        }

        /**
         * Returns the critical density (density of a flat universe) at redshift z in h^2 Msun/Mpc^3.
         */
        function rhoCritical(z){

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;
            let y   = Om0 * zp1 * zp1 * zp1 + Ode0 * Math.pow(zp1, 3 + 3*w) + Ok0 * zp1 * zp1;
            
            return RHO_CRIT0_ASTRO * y; 
        }

        /**
         * Returns the CMB temperature in Kelvin at redshift z.
         */
        function Tcmb(z){

            return Tcmb0 * (z + 1);
        }

        /**
         * Returns the dark-energy equation of state parameter at redshift z.
         */
        function wde(z){

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;
            
            return w;
        }

        /**
         * Returns the Hubble parameter at redshift z, in units of its present value.
         */
        function E(z){

            let zp1 = z + 1;
            let w   = w0 + wa * z / zp1;
            let y   = Om0 * zp1 * zp1 * zp1 + Ode0 * Math.pow(zp1, 3 + 3*w) + Ok0 * zp1 * zp1;
            
            return Math.sqrt(y);
        }

        /**
         * Returns the Hubble parameter at redshift z, in km/sec/Mpc.
         */
        function H(z){

            return h * 100.0 * E(z);
        }

        /**
         * Returns the logarithmic derivative of the Hubble parameter w.r.to the scale factor.
         */
        function dlnEdlnzp1(z){

            let zp1 = z + 1;
            let t;

            let x = Om0 * zp1 * zp1 * zp1;
            let y = 3*x;

            if (!flat){
                t  = Ok0 * zp1 * zp1;
                x += t;
                y += 2*t;
            }

            let w = w0 + wa * z / zp1;
            let f = 3*(1 + w) + 3 * wa / zp1 * Math.log(zp1);

            t  = Ode0 * Math.pow(zp1, 3 + 3*w);
            x += t;
            y += f*t;

            return 0.5 * y / x;
        }

        // integrals of z-functions

        /**
         * Returns the integral of (z+1)/E(z)^3 from redshift z_a to z_b.
         */
        function zIntegral_zp1_over_Ez3(za, zb){

            function logfunc(lnzp1){

                let zp1 = Math.exp(lnzp1);
                let w   = w0 + wa * (zp1 - 1) / zp1;
                let y   = Om0 * zp1 * zp1 * zp1 + Ode0 * Math.pow(zp1, 3 + 3*w) + Ok0 * zp1 * zp1;

                return zp1 * zp1 / Math.pow(y, 1.5);
            }

            let y = numeric.simpson(logfunc, Math.log(za + 1), Math.log(zb + 1), 1e-08, 100)
            return y;
        }

        /**
         * Returns the integral of 1/(z+1)E(z) from redshift z_a to z_b.
         */
        function zIntegral_1_over_zp1_Ez(za, zb){

            function logfunc(lnzp1){

                let zp1 = Math.exp(lnzp1);
                let w   = w0 + wa * (zp1 - 1) / zp1;
                let y   = Om0 * zp1 * zp1 * zp1 + Ode0 * Math.pow(zp1, 3 + 3*w) + Ok0 * zp1 * zp1;

                return 1 / Math.sqrt(y);
            }

            let y = numeric.simpson(logfunc, Math.log(za + 1), Math.log(zb + 1), 1e-08, 100)
            return y;
        }

        /**
         * Returns the integral of 1/E(z) from redshift z_a to z_b.
         */
        function zIntegral_1_over_Ez(za, zb){

            function logfunc(lnzp1){

                let zp1 = Math.exp(lnzp1);
                let w   = w0 + wa * (zp1 - 1) / zp1;
                let y   = Om0 * zp1 * zp1 * zp1 + Ode0 * Math.pow(zp1, 3 + 3*w) + Ok0 * zp1 * zp1;

                return zp1 / Math.sqrt(y);
            }

            let y = numeric.simpson(logfunc, Math.log(za + 1), Math.log(zb + 1), 1e-08, 100)
            return y;
        }

        // time

        /**
         * Returns the hubble time at redshift z in years.
         */
        function hubbleTime(z){

            let Hz = H(z) * (1e+12 / MPC * YEAR); // hubble paramtere in 1/yr 
            return 1 / Hz;
        }

        /**
         * Returns the age of the universe at redshift z in years.
         */
        function age(z){

            let inf = 1e+08;
            return zIntegral_1_over_zp1_Ez(z, inf) * hubbleTime(z);
        }

        /**
         * Returns the lookback time at redshift z in years. It is the difference between the present 
         * age of the universe and the age at z.
         */
        function lookbackTime(z){

            return age(0) - age(z);
        }

        // distances

        /**
         * Returns the comoving distance corresponding to the redshift z in Mpc.
         */
        function comovingDistance(z){

            let fac = C_SI / (100.0*h) / 1000.0;
            return fac * zIntegral_1_over_Ez(0.0, z);
        }

        /**
         * Returns the comoving coordinate distance corresponding to the redshift z in Mpc.
         */
        function comovingCoordinate(z){

            let x = comovingDistance(z);

            if (Math.abs(Ok0) < 1e-08){

                let k = Math.sqrt( Math.abs(Ok0) ) * ( 1e+05 * h / C_SI );

                if (Ok0 < 0){
                    return Math.sin(k*x) / k; // closed
                }
                return Math.sinh(k*x) / k; // open
            }
            return x; // flat
        }

        /**
         * Returns the angular diameter distance corresponding to the redshift z, in Mpc.
         */
        function angularDiameterDistance(z){

            let r = comovingCoordinate(z);
            return r / (z + 1);
        }

        /**
         * Returns the luminocity distance corresponding to the redshift z, in Mpc.
         */
        function luminocityDistance(z){

            let r = comovingCoordinate(z);
            return r * (z + 1);
        }

        // horizons

        /**
         * Returns the event horizon at redshift z in Mpc. Event horizon is the distance at which 
         * light emitted now can reach in future.
         */
        function eventHorizon(z){

            let c  = C_SI / 1000.0; // speed of light in km/sec
            let H0 = 100.0*h;

            let y = zIntegral_1_over_Ez(-0.99999999, z);
            return y * c / H0;
        }

        /**
         * Returns the particle horizon at redshift z in Mpc. Particle horizon is the distance from 
         * which the light detected now was emitted.
         */
        function particleHorizon(z){

            let c = C_SI / 1000.0; // speed of light in km/sec
            let H0 = 100.0*h;

            let inf = 1e+08;
            let y   = zIntegral_1_over_Ez(z, inf);
            return y * c / H0;
        }
        
        /**
         * Returns the Hubble horizon at redshift z in Mpc.
         */
        function hubbleHorizon(z){

            let c = C_SI / 1000.0; // speed of light in km/sec
            return c / H( z )      // Mpc
        }

        // linear growth factors

        /**
         * Tells if to use an exactly computed linear growth factor or its approximation by Carroll 
         * et al (1992).
         */
        function useExactGrowth(state){

            exactgf = state;
        }

        /**
         * Returns the linear growth factor of the density perturbations at redshift z, w.r.to that in 
         * a matter dominated universe.
         */
        function growth(z){

            function carroll92(z){

                let zp1 = z + 1;
                let w   = w0 + wa * z / zp1;

                let Omz, Odez;
                Omz   = Om0 * zp1 * zp1 * zp1;
                Odez  = Ode0 * Math.pow(zp1, 3 + 3*w);

                let y = Omz + Odez + Ok0 * zp1 * zp1;
                Omz   = Omz / y; 
                Odez  = Odez / y;

                let g = 2.5 * Omz / (Math.pow(Omz, 4.0/7.0) - Odez + ( 1 + Omz/2 ) * ( 1 + Odez/70 ));
                return g;
            }

            function growthExact(z){

                let y = zIntegral_zp1_over_Ez3(z, 1e+08);
                return 2.5 * Om0 * (z + 1) * E(z) * y;
            }

            if (exactgf){
                return growthExact(z);
            }
            return carroll92(z);
        }

        /**
         * Returns the linear growth factor of the density perturbations at redshift z.
         */
        function dplus(z){

            return growth(z) / (z + 1);
        }

        /**
         * Returns the logarithmic rate of linear growth factor w.r.to the scale factor, at redshift z.
         */
        function growthRate(z){

            if(exactgf){

                return 2.5 * Om(z) / growth(z, true) - dlnEdlnzp1(z);
            }
            return Math.pow(Om(z), 5.0/9.0);
        }

        // linear matter power spectrum

        /**
         * Set the linear matter power spectrum model. The value should be any of `eisenstein98_zb`, 
         * `eisenstein98_bao`, `eisenstein98_mdm` (not implemented yet) or `sugiyama95`. 
         */
        function setLinearPowerSpectrumModel(model){

            ps.setActive(model);
            createPkTable(kquadargs);
        }

        /**
         * Set the filter function model to use. Avaliable models are `tophat` and `gauss`. 
         */
        function setFilter(model){

            _filter = model;
        }

        /**
         * Returns the linear transfer function.
         * 
         * @param k Wavenumber in h/Mpc.
         * @param z Redshift.
         */
        function transfer(k, z){

            return ps.transfer(k, z);
        }

        /**
         * Returns the linear matter power spectrum, or if the `full` option is set, an object with 
         * fields `k` (wavenumber), `Pk` (linear power), `Dk` (dimensionless power) and `Tk` (transfer 
         * function).
         * 
         * @param k Wavenumber in h/Mpc.
         * @param z Redshift.
         * @param dim If true, return the dimensionless matter power spectrum.
         * @param full If true, returns an object containing the all power spectrum values.
         */
        function matterPowerSpectrum(k, z, dim, full){
            
            if (z == undefined){ z = 0.0; }
            if (dim == undefined){ dim = true; }
            if (full == undefined){ full = false; }

            let tk = transfer(k, z);
            let pk = A * Math.pow(k, ns) * tk * tk;

            var dz = 1.0;
            if (z > 0){
                
                dz = dplus(z) / dplus(0);
                pk = pk * dz * dz;
            }

            let dk = Math.pow(k, 3) * pk * ONE_2SQRPI;

            if (full){
                return { k: k, Pk: pk, Tk: tk, Dk: dk };
            }

            if (dim == false){
                return dk;
            }
            return pk;
        }

        /**
         * Returns the value of the 2-point correlation function.
         * 
         * @param r Distance or seperation between the points in Mpc/h.
         * @param z Redshift.
         */
        function matterCorrelation(r, z){

            let n    = Pktab.length;
            let lnka = Pktab[0][0], lnkb = Pktab[n-1][0];
            let dlnk = (lnkb - lnka) / (n-1);
            let ys   = [];
            for (let i = 0; i < n; i++){

                let k   = Math.exp(Pktab[i][0])
                let pk  = Pktab[i][1];
                
                let kr = k*r;
                let j0 = Math.sin(kr) / kr;

                // apply extra smoothing to remove small scale oscillation effectes
                let smooth = 0.01;
                smooth     = Math.exp( -kr * kr * smooth * smooth );

                ys.push( pk * j0 * smooth );
            }

            let y = numeric.simpsonSum(ys, dlnk);


            var dz = 1.0;
            if (z > 0){
                
                dz = dplus(z) / dplus(0);
                y  = y * dz * dz;
            }

            return A * y;
        }

        /**
         * Returns the smoothed linear matter variance of the density perturbations.
         * 
         * @param r Smoothing scale in Mpc/h.
         * @param z Redshift.
         */
        function variance(r, z){

            let n    = Pktab.length;
            let lnka = Pktab[0][0], lnkb = Pktab[n-1][0];
            let dlnk = (lnkb - lnka) / (n-1);
            let ys   = [];
            for (let i = 0; i < n; i++){

                let k   = Math.exp(Pktab[i][0])
                let pk  = Pktab[i][1];
                let wk  = filter(k*r, 0, _filter);

                ys.push( pk * wk * wk );
            }

            let y = numeric.simpsonSum(ys, dlnk);


            var dz = 1.0;
            if (z > 0){
                
                dz = dplus(z) / dplus(0);
                y  = y * dz * dz;
            }

            return A * y;
        }

        /**
         * Returns the logarithmic derivative of the variance w.r.to the smoothing radius.
         * 
         * @param r Smoothing scale in Mpc/h.
         * @param z Redshift.
         */
        function dlnsdlnr(r, z){

            let n    = Pktab.length;
            let lnka = Pktab[0][0], lnkb = Pktab[n-1][0];
            let dlnk = (lnkb - lnka) / (n-1);
            let ys   = [];
            for (let i = 0; i < n; i++){

                let k   = Math.exp(Pktab[i][0])
                let pk  = Pktab[i][1];
                let wk  = filter(k*r, 0, _filter);
                let dwk = filter(k*r, 1, _filter);

                ys.push( 2 * k * pk * wk * dwk );
            }

            let y = numeric.simpsonSum(ys, dlnk);

            var dz = 1.0;
            if (z > 0){
                
                dz = dplus(z) / dplus(0);

                y = y * dz * dz;
            }

            return 0.5 * r * A * y / variance(r, z);
        }

        /**
         * To do
         */
        function d2lnsdlnr2(){} // to do

        /**
         * Returns the effective slope of the power spectrum at wavenumber k.
         * 
         * @param k Wavenumber in h/Mpc.
         * @param z Redshift.
         * @param relstep Relative stepsize for differentiation. Default is 0.01. 
         */
        function effectiveIndex(k, z, relstep){

            if (relstep == undefined){
                relstep = 0.01;
            }

            let k1 = k * (1 - relstep), k2 = k * (1 + relstep);
            let p1 = matterPowerSpectrum(k1, z, true, false);
            let p2 = matterPowerSpectrum(k2, z, true, false);

            let dlnp = Math.log(p2) - Math.log(p1);
            let dlnk = Math.log(k2) - Math.log(k1);

            return dlnp / dlnk;
        }

        /**
         * Returns the smoothing radius corresponding to the variance value.
         * 
         * @param sigma Value of the variance.
         * @param z Redshift
         */
        function radius(sigma, z){

            // solve for radius - newton's method
            let v  = Math.pow( sigma / dz, 2 );

            let r0 = v;
            let not_converged = true;

            let v0, dv0, r;
            while (not_converged){

                v0 = variance(r0, z);

                let r1 = r0 * (1 - 0.01);
                let r2 = r0 * (1 + 0.01);
                let v1 = variance(r1, z);
                let v2 = variance(r2, z);

                dv0 = (v2 - v1) / (r2 - r1);
                
                r = r0 - v0 / dv0;

                if (Math.abs(1 - r0 / r) < 1e-05){
                    not_converged = false;
                }

                r0 = r;
            }

            return r0;
        }

        /** 
         * Normalize the linear matter power spectrum so that the sigma-8 parameter matches with the
         * calculated value.
         */
        function normalize(){

            A = 1.0;

            let var8 = variance(8.0, 0.0);
            A        = sigma8 * sigma8 / var8;
        }

        /**
         * Returns the normalisation of the matter power spectrum.
         */
        function powerSpectrumNorm(){

            return A;
        }

        // power spectrum table

        /**
         * Set integration settings for k-space integrations. The options should be given as an object 
         * with fields `ka` (lower limit in h/Mpc), `kb` (upper limit in h/Mpc) and `n` (the number of 
         * function evaluations). This also setup a power spectrum table for value lookup.
         */
        function setupPowerQuadArgs(args){

            kquadargs.ka = args.ka || kquadargs.ka;
            kquadargs.kb = args.kb || kquadargs.kb;
            kquadargs.n  = args.n  || kquadargs.n;

            return createPkTable(kquadargs);
        }

        /**
         * Create a power spectrum table to use for integration.
         */
        function createPkTable(args){

            let ka = args.ka;
            let kb = args.kb;
            let n  = args.n;

            let lnka = Math.log(ka), lnkb =Math.log(kb);
            let dlnk = (lnkb - lnka) / (n-1);
            let ys = [];
            for (let i = 0, lnk = lnka; i < n; i++, lnk += dlnk){

                let k  = Math.exp(lnk);
                let tk = transfer(k, 0);
                let Dk = Math.pow(k, ns + 3) * tk * tk / ( 2 * Math.PI * Math.PI ); // unnormalised dimensionless power  

                ys.push( [lnk, Dk] );
            }

            Pktab = ys;
        }

        // halo 

        /**
         * Returns the Lagrangian (comoving) radius in Mpc/h corresponding to the mass m (in Msun/h).
         */
         function lagrangianR(m){

            let r = Math.pow( 0.75*m / ( Math.PI * rho_m(0) ), 1/3 );
            return r;
        }

        /**
         * Returns the mass in Msun/h, corresponding to the Lagrangian radius r in Mpc/h.
         */
        function lagrangianM(r){

            let m = ( 4*Math.PI / 3.0 ) * Math.pow(r, 3) * rho_m(0);
            return m;
        }

        
        /** Public API for the cosmology module. */
        var publicAPI = {
                            parameters                 : model,
                            Om                         : Om,
                            Ob                         : Ob,
                            Onu                        : Onu,
                            Ode                        : Ode,
                            Ok                         : Ok,
                            rho_m                      : rho_m,
                            rho_de                     : rho_de,
                            rhoCritical                : rhoCritical,
                            Tcmb                       : Tcmb, // cmb temperature
                            wde                        : wde, // equation of state for dark-energy
                            H                          : H,
                            E                          : E,
                            zIntegral_1_over_Ez        : zIntegral_1_over_Ez,
                            zIntegral_1_over_zp1_Ez    : zIntegral_1_over_zp1_Ez,
                            zIntegral_zp1_over_Ez3     : zIntegral_zp1_over_Ez3, 
                            hubbleTime                 : hubbleTime,
                            lookbackTime               : lookbackTime,
                            age                        : age,
                            comovingDistance           : comovingDistance,
                            comovingCoordinate         : comovingCoordinate,
                            angularDiameterDistance    : angularDiameterDistance,
                            luminocityDistance         : luminocityDistance,
                            hubbleHorizon              : hubbleHorizon,
                            eventHorizon               : eventHorizon,
                            particleHorizon            : particleHorizon,
                            useExactGrowth             : useExactGrowth,
                            growth                     : growth,
                            growthRate                 : growthRate,
                            dplus                      : dplus,
                            transfer                   : transfer,
                            matterPowerSpectrum        : matterPowerSpectrum,
                            matterCorrelation          : matterCorrelation,
                            powerSpectrumNorm          : powerSpectrumNorm,
                            effectiveIndex             : effectiveIndex,
                            normalize                  : normalize,
                            variance                   : variance,
                            radius                     : radius,
                            dlnsdlnr                   : dlnsdlnr,
                            d2lnsdlnr2                 : d2lnsdlnr2,
                            lagrangianM                : lagrangianM,
                            lagrangianR                : lagrangianR,
                            massFunction               : massFunction,
                            linearBias                 : linearBias,
                            setupPowerQuadArgs         : setupPowerQuadArgs,
                            setFilter                  : setFilter,
                            setHaloMassfunctionModel   : setHaloMassfunctionModel,
                            setLinearPowerSpectrumModel: setLinearPowerSpectrumModel,
                            setLinearBiasModel         : setLinearBiasModel,
                        };


        // halo mass function

        var hmf = mass_function(publicAPI);
        hmf.setModel('tinker08');

        /**
         * Set the active halo mass-function model.
         */
        function setHaloMassfunctionModel(model){
            
            return hmf.setModel(model);
        }

        /**
         * Returns the halo mass-function in the format specified by the `fmt` argument, or an object 
         * containing the fields `m` (mass in Msun/h), `r` (Lagrangian radius in Mpc/h), `dndm`, `dndlnm` 
         * `f` (fitting function), `s` (variance/sigma) and `dlnsdlnm` (logarithmic derivative of the 
         * variance w.r.to mass).
         * 
         * @param m Mass of the halo in Msun/h.
         * @param z Redshift.
         * @param overdensity Overdensity value w.r.to mean background density. 
         * @param fmt Output format - any of `dndm`, `dndlnm`, `f` or `full`.
         */
        function massFunction(m, z, overdensity, fmt){

            return hmf.massFunction(m, z, overdensity, fmt);
        }

        var bias = linear_bias(publicAPI);
        bias.setModel('tinker10');

        /**
         * Set the active linear halo bias model.
         */
        function setLinearBiasModel(model){

            return bias.setModel(model);
        }

        /**
         * Returns the linear bias using the active model set.
         * @param nu Input argument, equal to $\sigma / \delta_c$
         * @param z  Redshift
         * @param overdensity Overdensity value w.r.to the mean background density.
         * @param correct Whether to use correction. Only for model `seljak04`. 
         */
        function linearBias(nu, z, overdensity, correct){

            return bias.bias(nu, z, overdensity, correct);
        }

        /**
         * To do
         */
        function haloPowerSpectrum(){}

        /**
         * To do
         */
        function haloCorrelation(){}


        return publicAPI;

    }

    /** Predefined cosmologies. */
    const predefined = {

        /** Cosmology model with parameters from Plank (2015) survey. */
        plank15  : cosmology({ h: 0.6790, Om0: 0.3065, Ob0: 0.0483, Ode0: 0.6935, sigma8: 0.8154, ns: 0.9681, Tcmb0: 2.7255, flat: false }),

        /** Cosmology model with parameters from Plank (2018) survey. */
        plank18: cosmology({ h: 0.6736, Om0: 0.3153, Ob0: 0.0493, Ode0: 0.6947, sigma8: 0.8111, ns: 0.9649, Tcmb0: 2.7255, flat: false }),

        /** Cosmology model with parameters from WMAP (2008) survey. */
        wmap08: cosmology({ h: 0.719, Om0: 0.2581, Ob0: 0.0441, Ode0: 0.7420, sigma8: 0.796, ns: 0.963, Tcmb0: 2.7255, flat: false }),

        /** Cosmology model for millanium simulation. */
        millanium: cosmology({ h: 0.73, Om0: 0.25, Ob0: 0.045, sigma8: 0.9, ns: 1.0, Tcmb0: 2.7255, flat: true }),

    };


    /** cosmology public API */
    return {
                cosmology     : cosmology,
                predefined    : predefined,
                power_spectrum: power_spectrum,
                mass_function : mass_function,
                linear_bias   : linear_bias,
           }
}

const cosmo = cosmologyModule();
