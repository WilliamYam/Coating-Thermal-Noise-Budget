function iopFilter = iopdownsamplingfilters(freq,sampleFreqString,filterCalcMethod)

% Coefficients copied from lines 289 to 354 of RCG 2.6.2 tagged release of
% <https://redoubt.ligo-wa.caltech.edu/svn/advLigoRTS/tags/advLigoRTS-2.6.2/src/fe/controller.c>
%
% Parsed using method outlined in 
% <https://gist.github.com/tobin/2206443>

% // All systems not running at 64K require up/down sampling to communicate I/O data
% // with IOP, which is running at 64K.
% // Following defines the filter coeffs for these up/down filters.
% #ifdef OVERSAMPLE

% filterCalcMethod = 'biquad';
% % filterCalcMethos = 'df2';
% sampleFreqString = '4k';

switch filterCalcMethod
    case 'biquad'
        % /* Recalculated filters in biquad form */
        % /* Oversamping base rate is 64K */
        
        switch sampleFreqString
            case '32k'
                % /* Coeffs for the 2x downsampling (32K system) filter */
                % static double __attribute__ ((unused)) feCoeff2x[9] =
                %        {0.053628649721183,
                %         0.2568759660371100, -0.3225906481359000, 1.2568801238621801, 1.6774135096891700,
                %        -0.2061764045745400, -1.0941543149527400, 2.0846376586498803, 2.1966597482716801};
                gain = 0.053628649721183;
                sos  = [ 0.2568759660371100, -0.3225906481359000, 1.2568801238621801, 1.6774135096891700;...
                        -0.2061764045745400, -1.0941543149527400, 2.0846376586498803, 2.1966597482716801];
                
            case '16k'
                % /* Coeffs for the 4x downsampling (16K system) filter */
                % static double __attribute__ ((unused)) feCoeff4x[9] =
                %        {0.014805052402446,
                %         0.7166258547451800, -0.0683289874517300, 0.3031629575762000, 0.5171469569032900,
                %         0.6838596423885499, -0.2534855521841101, 1.6838609161411500, 1.7447155374502499};
                gain = 0.014805052402446;
                sos  = [ 0.7166258547451800, -0.0683289874517300, 0.3031629575762000, 0.5171469569032900;...
                         0.6838596423885499, -0.2534855521841101, 1.6838609161411500, 1.7447155374502499];
                
            case '4k'
                % // New Brian Lantz 4k decimation filter
                % static double __attribute__ ((unused)) feCoeff16x[9] =
                %       {0.010203728365,
                %        0.8052941009065100, -0.0241751519071000, 0.3920490703701900, 0.5612099784288400,
                %        0.8339678987936501, -0.0376022631287799, -0.0131581721533700, 0.1145865116421301};
                gain = 0.010203728365;
                sos  = [ 0.8052941009065100, -0.0241751519071000, 0.3920490703701900, 0.5612099784288400;...
                         0.8339678987936501, -0.0376022631287799, -0.0131581721533700, 0.1145865116421301];
                
            case '2k'
                % /* Coeffs for the 32x downsampling filter (2K system) */
                % /* Original Rana coeffs from 40m lab elog */
                % static double __attribute__ ((unused)) feCoeff32x[9] =
                %        {0.0001104130574447,
                %         0.9701834961388200, -0.0010837026165800, -0.0200761119821899, 0.0085463156103800,
                %         0.9871502388637901, -0.0039246182095299, 3.9871502388637898, 3.9960753817904697};
                gain = 0.0001104130574447;
                sos  = [ 0.9701834961388200, -0.0010837026165800, -0.0200761119821899, 0.0085463156103800;...
                         0.9871502388637901, -0.0039246182095299, 3.9871502388637898, 3.9960753817904697];
        end
        
        % Put it into standard form
        
        % BiQuad coefficients are reported in rows of [a11, a12, c1, and c2]
        % instead of rows of [a1, a2, b1, b2] for DF2. So, we need to
        % compute a1, a2, b1, and b2. The relation is derived from pg 11 of
        % G0900928-v1
        rows = size(sos,1);
        
        b0 =  ones(rows,1);
        b1 = -sos(:,1) + sos(:,3) - 1;
        b2 =  sos(:,1) - sos(:,2) - sos(:,3) + sos(:,4);
        a0 =  ones(rows,1);
        a1 = -sos(:,1) - 1;
        a2 =  sos(:,1) - sos(:,2);        
        
    case 'df2'
        switch sampleFreqString
            % /* Oversamping base rate is 64K */
            case '32k'
                % /* Coeffs for the 2x downsampling (32K system) filter */
                % static double __attribute__ ((unused)) feCoeff2x[9] =
                %        {0.053628649721183,
                %        -1.25687596603711,    0.57946661417301,    0.00000415782507,    1.00000000000000,
                %        -0.79382359542546,    0.88797791037820,    1.29081406322442,    1.00000000000000};
                gain = 0.053628649721183;
                sos  = [-1.25687596603711,    0.57946661417301,    0.00000415782507,    1.00000000000000;...
                        -0.79382359542546,    0.88797791037820,    1.29081406322442,    1.00000000000000];
                
            case '16k'
                % /* Coeffs for the 4x downsampling (16K system) filter */
                % static double __attribute__ ((unused)) feCoeff4x[9] =
                %        {0.014805052402446,
                %        -1.71662585474518,    0.78495484219691,   -1.41346289716898,   0.99893884152400,
                %        -1.68385964238855,    0.93734519457266,    0.00000127375260,   0.99819981588176};
                gain = 0.014805052402446;
                sos  = [-1.71662585474518,    0.78495484219691,   -1.41346289716898,   0.99893884152400;...
                        -1.68385964238855,    0.93734519457266,    0.00000127375260,   0.99819981588176];
        
            case '4k'
                % // New Brian Lantz 4k decimation filter
                % static double __attribute__ ((unused)) feCoeff16x[9] =
                %        {0.010203728365,
                %        -1.80529410090651,   0.82946925281361,  -1.41324503053632,   0.99863016087226,
                %        -1.83396789879365,   0.87157016192243,  -1.84712607094702,   0.99931484571793};
                gain = 0.010203728365;
                sos = [-1.80529410090651,   0.82946925281361,  -1.41324503053632,   0.99863016087226;...
                       -1.83396789879365,   0.87157016192243,  -1.84712607094702,   0.99931484571793];
                

%     #if 0
%         % /* Coeffs for the 32x downsampling filter (2K system) */
%         % /* Original Rana coeffs from 40m lab elog */
%         % static double feCoeff32x[9] =
%                 {0.0001104130574447,
%                 -1.97018349613882,    0.97126719875540,   -1.99025960812101,    0.99988962634797,
%                 -1.98715023886379,    0.99107485707332,    2.00000000000000,    1.00000000000000};
%     #endif
            case '2k'
                % /* Coeffs for the 32x downsampling filter (2K system) per Brian Lantz May 5, 2009 */
                % static double __attribute__ ((unused)) feCoeff32x[9] =
                %        {0.010581064947739,
                %        -1.90444302586137,    0.91078434629894,   -1.96090276933603,    0.99931924465090,
                %        -1.92390910024681,    0.93366146580083,   -1.84652529182276,    0.99866506867980};
                gain = 0.010581064947739;
                sos = [-1.90444302586137,    0.91078434629894,   -1.96090276933603,    0.99931924465090;...
                       -1.92390910024681,    0.93366146580083,   -1.84652529182276,    0.99866506867980];
        end
        
        % Put it into standard form
        %
        % The coefficients stored in the C code assume a form of the second order
        % section where the leading coefficient in the numerator and the
        % denominator are both 1. Here we insert those two 1's into the matrix,
        % and also swap the numerator and denominator.
        rows = size(sos,1);
        
        b0 = ones(rows,1);
        b1 = sos(:,3);
        b2 = sos(:,4);
        a0 = ones(rows,1);
        a1 = sos(:,1);
        a2 = sos(:,2); 
end

iopFilter.sosCoefficients = [ b0 b1 b2 a0 a1 a2 ];

omega = 2*pi*freq;

iopFilter.f  = gain * sos2freqresp(iopFilter.sosCoefficients, omega, 2^16);
[A,B,C,D] = sos2ss(iopFilter.sosCoefficients,gain);
iopFilter.ss = d2c(ss(A,B,C,D,1/2^16), 'tustin');