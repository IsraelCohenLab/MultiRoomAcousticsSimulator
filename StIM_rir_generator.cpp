/*
Program     : StIM Room Impulse Response Generator
 
Description : Computes the response of an acoustic source to a microphone
              residing in an adjacent room using the Structuralimage method [1,2].
 
              [1] J.B. Allen and D.A. Berkley,
              Image method for efficiently simulating small-room acoustics,
              Journal Acoustic Society of America, 65(4), April 1979, p 943.
 
              [2] E.Shalev I.Cohen and D.Levov,
              ndoors audio classification with structure image methodfor simulating multi-room acoustics'
              The Journal of the Acoustical Society of America, 150(4):3059–3073, 2021.
*/

#define _USE_MATH_DEFINES

#include "matrix.h"
#include "mex.h"
#include "math.h"

#define ROUND(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

#ifndef M_PI 
    #define M_PI 3.14159265358979323846 
#endif


/*
    This function receives:
        L:          Dimention of a cubic box\room. The origins is assumed to be located in the (0,0,0) location.
        xv:         The (x,y,z) coordination of the virtual source\receiver. 
        xf:         The (x,y,z) coordination of the fixed source\receiver. 
    The function returns:
        face_num:   The numbe of the face through which the ray penetrates. 0 represents no penetration.
                    1 and 2 represent the wall in x=0 and x=L_x respectively.
                    3 and 4 represent the wall in y=0 and y=L_y respectively.
                    5 and 6 represent the wall in z=0 and z=L_z respectively.
*/

int box_ray(const double L[],double xv[],const double xf[]){
    
    // Negated deltas
    double ndx = xv[0] - xf[0];
    double ndy = xv[1] - xf[1];
    double ndz = xv[2] - xf[2];
    
    // Sizes scaled by the negated deltas
    double sxy = ndx * L[1];
    double sxz = ndx * L[2];
    double syx = ndy * L[0];
    double syz = ndy * L[2];
    double szx = ndz * L[0];
    double szy = ndz * L[1];
    // Cross terms
    double cxy = xf[0]*xv[1] - xf[1]*xv[0];
    double cxz = xf[0]*xv[2] - xf[2]*xv[0];
    double cyz = xf[1]*xv[2] - xf[1]*xv[2];

    // Absolute delta products
    double axy = abs(ndx*ndy);
    double axz = abs(ndx*ndz);
    double ayz = abs(ndy*ndz);
    double axyz = abs(ndz*axy);
    

    // Default to "no intersection"
    int face_num = 0;
    float face_tau = abs(ndz*axy);
    

    
    // These variables are no longer used:
    double tau;
    
    if (xv[0] < 0 and 0 < xf[0]){
        // Face 1: x == 0
        tau = -xv[0] * ayz;
        if (tau < face_tau && cxy >= 0 && cxz >= 0 && cxy <= -sxy && cxz <= -sxz){
            face_tau = tau;
            face_num = 1;
        }
    }
    
    else if (xf[0] < L[0] && L[0] < xv[0]){
        // Face 2: x == size[0]
        tau = (xv[0] - L[0]) * ayz;
        if (tau < face_tau && cxy <= syx && cxz <= szx && cxy >= syx - sxy && cxz >= szx - sxz){
            face_tau = tau;
            face_num = 2;
        }
    }

    if (xv[1] < 0 && xf[1] > 0){
        // Face 3: y == 0
        tau = -xv[1] * axz;
        if (tau < face_tau && cxy <= 0 && cyz >= 0 && cxy >= syx && cyz <= -syz){
            face_tau = tau;
            face_num = 3;
        }
    }

    else if (xv[1] > L[1] && xf[1] < L[1]){
        // Face 4: y == size[1]
        tau = (xv[1] - L[1]) * axz;
        if (tau < face_tau && cxy >= -sxy && cyz <= szy && cxy <= syx - sxy && cyz >= szy - syz){
            face_tau = tau;
            face_num = 4;
        }
    }

    if (xv[2] < 0 && xf[2] > 0){
        // Face 5: z == 0
        tau = -xv[2] * axy;
        if (tau < face_tau && cxz <= 0 && cyz <= 0 && cxz >= szx && cyz >= szy){
            face_tau = tau;
            face_num = 5;
        }
    }

    else if (xv[2] > L[2] && xf[2] < L[2]){
        // Face 6: z == size[2]
        tau = (xv[2] - L[2]) * axy;
        if (tau < face_tau && cxz >= -sxz && cyz >= -syz && cxz <= szx - sxz && cyz <= szy - syz){
            face_tau = tau;
            face_num = 6;
        }
    }
    
    return face_num;
}


double sinc(double x)
{
    if (x == 0)
        return(1.);
    else
        return(sin(x)/x);
}


double sim_microphone(double x, double y, double z, double* angle, char mtype)
{
    if (mtype=='b' || mtype=='c' || mtype=='s' || mtype=='h')
    {
        double gain, vartheta, varphi, rho;

        // Polar Pattern         rho
        // ---------------------------
        // Bidirectional         0
        // Hypercardioid         0.25    
        // Cardioid              0.5
        // Subcardioid           0.75
        // Omnidirectional       1

        switch(mtype)
        {
        case 'b':
            rho = 0;
            break;
        case 'h':
            rho = 0.25; 
            break;
        case 'c':
            rho = 0.5;
            break;
        case 's':
            rho = 0.75;
            break;
        };
                
        vartheta = acos(z/sqrt(pow(x,2)+pow(y,2)+pow(z,2)));
        varphi = atan2(y,x);

        gain = sin(M_PI/2-angle[1]) * sin(vartheta) * cos(angle[0]-varphi) + cos(M_PI/2-angle[1]) * cos(vartheta);
        gain = rho + (1-rho) * gain;
                
        return gain;
    }
    else
    {
        return 1;
    }
}


/*  
    Reseiver image method function computes the room impulse response of a fixed source, outside a room. 
    The functin uses receiver imageing and not source imaging, as in the original image method.
    If the original receiver's room, is not on the line between a virtual receiver and the fixed source, the path is omitted as an artifect. 
 */
void ReciverIM(double c, double fs,const double* rr, int nMicrophones, const double* ss, const double* LL ,const double* beta_input, int nDimension, int nOrder, char* microphone_type,
    int nSamples,double* angle, double* imp)
{
    
    double*         beta = new double[6];
    double          reverberation_time = 0;
    bool            flag;
    
    // Setting some variables:
    
    // Reflection coefficients or reverberation time?
    if (sizeof(beta_input)/sizeof(beta_input[0])==1)
    {
        double V = LL[0]*LL[1]*LL[2];
        double S = 2*(LL[0]*LL[2]+LL[1]*LL[2]+LL[0]*LL[1]);
        reverberation_time = beta_input[0];
        if (reverberation_time != 0) {
            double alfa = 24*V*log(10.0)/(c*S*reverberation_time);
            if (alfa > 1)
                mexErrMsgTxt("Error: The reflection coefficients cannot be calculated using the current "
                             "room parameters, i.e. room size and reverberation time.\n           Please "
                             "specify the reflection coefficients or change the room parameters.");
            for (int i=0;i<6;i++)
                beta[i] = sqrt(1-alfa);
        } 
        else 
        {
            for (int i=0;i<6;i++)
                beta[i] = 0;
        }
    }
    else
    {
        for (int i=0;i<6;i++)
            beta[i] = beta_input[i];
    }
    
    if (nDimension == 2)
    {
     	beta[4] = 0;
     	beta[5] = 0;
    }
    
    // Temporary variables and constants (image-method)
    const double Fc = 1; // The cut-off frequency equals fs/2 - Fc is the normalized cut-off frequency.
    const int    Tw = 2 * ROUND(0.004*fs); // The width of the low-pass FIR equals 8 ms
    const double cTs = c/fs;
    double*      LPI = new double[Tw];
    double*      r = new double[3];
    double*      xp = new double[3];
    double*      s = new double[3];
    double*      L = new double[3];
    double       Rm[3];
    double       Rp_plus_Rm[3]; 
    double       refl[3];
    double       fdist,dist;
    double       gain;
    int          startPosition;
    int          n1, n2, n3;
    int          q, j, k;
    int          mx, my, mz;
    int          n;
    int          order;
    int          face_num;
    
    
    // Set spatial ampled variable
    s[0] = ss[0]/cTs; s[1] = ss[1]/cTs; s[2] = ss[2]/cTs;
    L[0] = LL[0]/cTs; L[1] = LL[1]/cTs; L[2] = LL[2]/cTs;
    
    // Loop all microphons in a microphone array
    for (int idxMicrophone = 0; idxMicrophone < nMicrophones ; idxMicrophone++)
    {
        
        // [x_1 x_2 ... x_N y_1 y_2 ... y_N z_1 z_2 ... z_N]
        r[0] = rr[idxMicrophone + 0*nMicrophones] / cTs;
        r[1] = rr[idxMicrophone + 1*nMicrophones] / cTs;
        r[2] = rr[idxMicrophone + 2*nMicrophones] / cTs;

        n1 = (int) ceil(nSamples/(2*L[0]));
        n2 = (int) ceil(nSamples/(2*L[1]));
        n3 = (int) ceil(nSamples/(2*L[2]));
        

        // Generate room impulse response
        // Loop all possible duplicagtions on the x,y, and z axes:
        for (mx = -n1 ; mx <= n1 ; mx++)// x axis
        {
            Rm[0] = 2*mx*L[0];

            for (my = -n2 ; my <= n2 ; my++)// y axis
            {
                Rm[1] = 2*my*L[1];

                for (mz = -n3 ; mz <= n3 ; mz++)// z axis
                {
                    Rm[2] = 2*mz*L[2];

                    // Loop all possible imaging variation for the axes:
                    for (q = 0 ; q <= 1 ; q++)//x axis
                    {
                        Rp_plus_Rm[0] = (1-2*q)*r[0] - s[0] + Rm[0];// x axis distance between virtual reciver and fixed source
                        xp[0] = 2*mx*LL[0] + (1-2*q)*rr[idxMicrophone + 0*nMicrophones];// virtual reciver x position
                        refl[0] = pow(beta[0], abs(mx-q)) * pow(beta[1], abs(mx));

                        for (j = 0 ; j <= 1 ; j++)// y axis
                        {
                            Rp_plus_Rm[1] = (1-2*j)*r[1] - s[1] + Rm[1];// y axis distance between virtual reciver and fixed source
                            xp[1] = 2*my*LL[1] + (1-2*j)*rr[idxMicrophone + 1*nMicrophones];// virtual reciver y position
                            refl[1] = pow(beta[2], abs(my-j)) * pow(beta[3], abs(my));

                            for (k = 0 ; k <= 1 ; k++)// z axis
                            {
                                Rp_plus_Rm[2] = (1-2*k)*r[2] - s[2] + Rm[2];// z axis distance between virtual reciver and fixed source
                                xp[2] = 2*mz*LL[2] + (1-2*k)*rr[idxMicrophone + 2*nMicrophones];// virtual reciver z position
                                refl[2] = pow(beta[4],abs(mz-k)) * pow(beta[5], abs(mz));
                                
                                order = abs(2*mx-q)+abs(2*my-j)+abs(2*mz-k);

                                dist = sqrt(pow(Rp_plus_Rm[0], 2) + pow(Rp_plus_Rm[1], 2) + pow(Rp_plus_Rm[2], 2)); // Audio traveling distance
                                
                                face_num = box_ray(LL, xp, ss); // Determine the penetration face number                        
                                
                                if (face_num == 0){
                                    continue;}// Do not proceed for artifact pathes
                                if (abs(2*mx-q)+abs(2*my-j)+abs(2*mz-k) <= nOrder || nOrder == -1)
                                {
                                    fdist = floor(dist);
                                    if (fdist < nSamples)
                                    {
                                        // Calculate gain with respect to traveling distance and reflection order
                                        gain = sim_microphone(Rp_plus_Rm[0], Rp_plus_Rm[1], Rp_plus_Rm[2], angle, microphone_type[0])
                                            * refl[0]*refl[1]*refl[2]/(4*M_PI*dist*cTs);
                                        
                                        // Set the correct gain for the penetration face
                                        gain = gain*(1- beta[face_num-1])/ beta[face_num-1];
                                        
                                        // Add to the total room impulse risponse
                                        for (n = 0 ; n < Tw ; n++)
                                            LPI[n] =  0.5 * (1 - cos(2*M_PI*((n+1-(dist-fdist))/Tw))) * Fc * sinc(M_PI*Fc*(n+1-(dist-fdist)-(Tw/2)));

                                        startPosition = (int) fdist-(Tw/2)+1;
                                        for (n = 0 ; n < Tw; n++){
                                            if (startPosition+n >= 0 && startPosition+n < nSamples){
                                                imp[idxMicrophone + nMicrophones*(startPosition+n)] += gain * LPI[n];
                                            }
                                        }// n for loop
                                    }// conditioning for number of samples
                                }// conditioning for none-artifact sources
                            }// k for loop
                        }// j for loop
                    }// q for loop
                }// mz for loop
            }// my for loop
        }// mx for loop
    }// Microphone for loop
    
    
    delete[] beta;
    delete[] LPI;            
    delete[] r;
    delete[] s;
    delete[] L;
    
}

/*  
    Main function called from MATLAB. Original image method function computes the room impulse response of a fixed reciever, outside a room. 
    If the original source's room, is not on the line between a virtual source and the fixed reciver, the path is omitted as an artifect. 
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Receive variables from MATLAB call
    if (nrhs == 0)
    {
        mexPrintf("--------------------------------------------------------------------\n"
            "| Room Impulse Response Generator                                  |\n"
            "|                                                                  |\n"
            "| Computes the response of an acoustic source to one or more       |\n"
            "| microphones in a reverberant room using the image method [1,2].  |\n"
            "|                                                                  |\n"
            "| [1] J.B. Allen and D.A. Berkley,                                 |\n"
            "|     Image method for efficiently simulating small-room acoustics,|\n"
            "|     Journal Acoustic Society of America,                         |\n"
            "|     65(4), April 1979, p 943.                                    |\n"
            "|                                                                  |\n"
            "| [2] E.Shalev I.Cohen and D.Levov,                                |\n"
            "|     Indoors audio classification with structure image method for |\n"
            "|     simulating multi-room acoustics,                             |\n"
            "|     The Journal of the Acoustical Society of America,            |\n"
            "|     150(4):3059–3073, 2021.                                      |\n"
            "--------------------------------------------------------------------\n\n"
            "function [h, beta_hat] = rir_generator(c, fs, r, s, Lr, Ls, beta_r, beta_s, nsample,\n"
            " mtype, order, dim, orientation, hp_filter);\n\n"
            "Input parameters:\n"
            " c           : sound velocity in m/s.\n"
            " fs          : sampling frequency in Hz.\n"
            " r           : 1 x 3 array specifying the (x,y,z) coordinates of the\n"
            "               receiver(s) in m.\n"
            " s           : 1 x 3 vector specifying the (x,y,z) coordinates of the\n"
            "               source in m.\n"
            " Lr          : 1 x 3 vector specifying the receiver's room dimensions (x,y,z) in m.\n"
            " Ls          : 1 x 3 vector specifying the source's room dimensions (x,y,z) in m.\n"
            " beta_r      : 1 x 6 vector specifying the reflection coefficients in the receiver's room\n"
            "               [beta_x1 beta_x2 beta_y1 beta_y2 beta_z1 beta_z2] or\n"
            "               beta = reverberation time (T_60) in seconds.\n"
            " beta_s      : 1 x 6 vector specifying the reflection coefficients in the source's room\n"
            "               [beta_x1 beta_x2 beta_y1 beta_y2 beta_z1 beta_z2] or\n"
            "               beta = reverberation time (T_60) in seconds.\n"
            " nsample     : number of samples to calculate, default is T_60*fs.\n"
            " mtype       : [omnidirectional, subcardioid, cardioid, hypercardioid,\n"
            "               bidirectional], default is omnidirectional.\n"
            " order       : reflection order, default is -1, i.e. maximum order.\n"
            " dim         : room dimension (2 or 3), default is 3.\n"
            " orientation : direction in which the microphones are pointed, specified using\n"
            "               azimuth and elevation angles (in radians), default is [0 0].\n"
            " hp_filter   : use 'false' to disable high-pass filter, the high-pass filter\n"
            "               is enabled by default.\n\n"
            "Output parameters:\n"
            " h           : M x nsample matrix containing the calculated room impulse\n"
            "               response(s).\n"
            " beta_hat    : In case a reverberation time is specified as an input parameter\n"
            "               the corresponding reflection coefficient is returned.\n\n");
    }
    else
    {
         mexPrintf("Cross-Room Impulse Response Generator (Version 1.1.2021) by Erez Shalev\n");
    }

    // Check for proper number of arguments
    if (nrhs < 8)
        mexErrMsgTxt("Error: There are at least eight input parameters required.");
    if (nrhs > 14)
        mexErrMsgTxt("Error: Too many input arguments.");
    if (nlhs > 2)
        mexErrMsgTxt("Error: Too many output arguments.");
    
    // Check for proper arguments
    if (!(mxGetN(prhs[0])==1) || !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Invalid input c arguments!");// C
    if (!(mxGetN(prhs[1])==1) || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
        mexErrMsgTxt("Invalid input fs arguments!");// fs
    if (!(mxGetN(prhs[2])==3) || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
        mexErrMsgTxt("Invalid input s arguments!");// s
    if (!(mxGetN(prhs[3])==3) || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]))
        mexErrMsgTxt("Invalid input r arguments!");// r
    if (!(mxGetN(prhs[4])==3) || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
        mexErrMsgTxt("Invalid input Lr arguments!");// Lr
    if (!(mxGetN(prhs[5])==3) || !mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]))
        mexErrMsgTxt("Invalid input Ls arguments!");// Ls
    if (!(mxGetN(prhs[6])==6 || mxGetN(prhs[6])==1) || !mxIsDouble(prhs[6]) || mxIsComplex(prhs[6]))
        mexErrMsgTxt("Invalid input beta_r arguments!");// beta_r
    if (!(mxGetN(prhs[7])==6 || mxGetN(prhs[7])==1) || !mxIsDouble(prhs[7]) || mxIsComplex(prhs[7]))
        mexErrMsgTxt("Invalid input beta_s arguments!");// beta_s

    // Load parameters
    double          c = mxGetScalar(prhs[0]);
    double          fs = mxGetScalar(prhs[1]);
    const double*   rr = mxGetPr(prhs[2]);
    int             nMicrophones = (int) mxGetM(prhs[2]);
    const double*   ss = mxGetPr(prhs[3]);
    const double*   LL_r = mxGetPr(prhs[4]);
    const double*   LL_s = mxGetPr(prhs[5]);
    const double*   beta_input_r = mxGetPr(prhs[6]);
    const double*   beta_input_s = mxGetPr(prhs[7]);
    double*         beta_s = new double[6];
    double*         beta = new double[6];
    int             nSamples;
    char*           microphone_type;
    int             nOrder;
    int             nDimension;
    double          angle[2];
    int             isHighPassFilter;
    double          reverberation_time = 0;
    bool            flag;

    // Reflection coefficients or reverberation time?
    if (mxGetN(prhs[7])==1)
    {
        double V = LL_s[0]*LL_s[1]*LL_s[2];
        double S = 2*(LL_s[0]*LL_s[2]+LL_s[1]*LL_s[2]+LL_s[0]*LL_s[1]);
        reverberation_time = beta_input_s[0];
        if (reverberation_time != 0) {
            double alfa = 24*V*log(10.0)/(c*S*reverberation_time);
            if (alfa > 1)
                mexErrMsgTxt("Error: The reflection coefficients cannot be calculated using the current "
                             "room parameters, i.e. room size and reverberation time.\n           Please "
                             "specify the reflection coefficients or change the room parameters.");
            for (int i=0;i<6;i++)
                beta[i] = sqrt(1-alfa);
        } 
        else 
        {
            for (int i=0;i<6;i++)
                beta[i] = 0;
        }
    }
    else
    {
        for (int i=0;i<6;i++)
            beta[i] = beta_input_s[i];
    }
    

    // High-pass filter (optional)
    if (nrhs > 13 &&  mxIsEmpty(prhs[13]) == false)
    {
        isHighPassFilter = (int) mxGetScalar(prhs[13]);
    }
    else
    {
        isHighPassFilter = 1;
    }

    // 3D Microphone orientation (optional)
    if (nrhs > 12 &&  mxIsEmpty(prhs[12]) == false)
    {
        const double* orientation = mxGetPr(prhs[12]);
        if (mxGetN(prhs[12]) == 1)
        {     
            angle[0] = orientation[0];
            angle[1] = 0;
        }
        else
        {
            angle[0] = orientation[0];
            angle[1] = orientation[1];
        }
    }
    else
    {
        angle[0] = 0;
        angle[1] = 0;
    }
    
    // Room Dimension (optional)
    if (nrhs > 11 &&  mxIsEmpty(prhs[11]) == false)
    {
        nDimension = (int) mxGetScalar(prhs[11]);
        if (nDimension != 2 && nDimension != 3)
            mexErrMsgTxt("Invalid input arguments!");

        if (nDimension == 2)
        {
            beta[4] = 0;
            beta[5] = 0;
        }
    }
    else
    {
        nDimension = 3;
    }

    // Reflection order (optional)
    if (nrhs > 10 &&  mxIsEmpty(prhs[10]) == false)
    {
        nOrder = (int) mxGetScalar(prhs[10]);
        if (nOrder < -1)
            mexErrMsgTxt("Invalid input arguments!");
    }
    else
    {
        nOrder = -1;
    }

    // Type of microphone (optional)
    if (nrhs > 9 &&  mxIsEmpty(prhs[9]) == false)
    {
        microphone_type = new char[mxGetN(prhs[9])+1];
        mxGetString(prhs[9], microphone_type, mxGetN(prhs[9])+1);
    }
    else
    {
        microphone_type = new char[1];
        microphone_type[0] = 'o';
    }

    // Number of samples (optional)
    if (nrhs > 8 &&  mxIsEmpty(prhs[8]) == false)
    {
        nSamples = (int) mxGetScalar(prhs[8]);
    }
    else
    {
        if (mxGetN(prhs[7])>1)
        {
            double V = LL_s[0]*LL_s[1]*LL_s[2];
            double alpha = ((1-pow(beta[0],2))+(1-pow(beta[1],2)))*LL_s[1]*LL_s[2] +
                ((1-pow(beta[2],2))+(1-pow(beta[3],2)))*LL_s[0]*LL_s[2] +
                ((1-pow(beta[4],2))+(1-pow(beta[5],2)))*LL_s[0]*LL_s[1];
            reverberation_time = 24*log(10.0)*V/(c*alpha);
            if (reverberation_time < 0.128)
                reverberation_time = 0.128;
        }
        nSamples = (int) (reverberation_time * fs);
    }

    // Create output vector
    plhs[0] = mxCreateDoubleMatrix(nMicrophones, nSamples, mxREAL);
    double* imp = mxGetPr(plhs[0]);

    // Temporary variables and constants (high-pass filter)
    const double W = 2*M_PI*100/fs; // The cut-off frequency equals 100 Hz
    const double R1 = exp(-W);
    const double B1 = 2*R1*cos(W);
    const double B2 = -R1 * R1;
    const double A1 = -(1+R1);
    double       X0;
    double*      Y = new double[3];

    // Temporary variables and constants (image-method)
    const double Fc = 1; // The cut-off frequency equals fs/2 - Fc is the normalized cut-off frequency.
    const int    Tw = 2 * ROUND(0.004*fs); // The width of the low-pass FIR equals 8 ms
    const double cTs = c/fs;
    double*      LPI = new double[Tw];
    double*      r = new double[3];
    double*      r_rIM = new double[3];// receiver location for receiver IM
    double*      xp = new double[3];
    double*      s = new double[3];
    double*      s_rIM = new double[3];// source location for receiver IM
    double*      L = new double[3];
    double       Rm[3];
    double       Rp_plus_Rm[3]; 
    double       refl[3];
    double       fdist,dist;
    double       gain;
    int          startPosition;
    int          n1, n2, n3;
    int          q, j, k;
    int          mx, my, mz;
    int          n;
    int          order;
    int          outsideSourceFlag;
    int          face_num = 0;
    
    // Set spatial ampled variable
    s[0] = ss[0]/cTs; s[1] = ss[1]/cTs; s[2] = ss[2]/cTs;
    L[0] = LL_s[0]/cTs; L[1] = LL_s[1]/cTs; L[2] = LL_s[2]/cTs;
    
    // Loop all microphons in a microphone array
    for (int idxMicrophone = 0; idxMicrophone < nMicrophones ; idxMicrophone++)
    {
        // [x_1 x_2 ... x_N y_1 y_2 ... y_N z_1 z_2 ... z_N]
        r[0] = rr[idxMicrophone + 0*nMicrophones] / cTs;
        r[1] = rr[idxMicrophone + 1*nMicrophones] / cTs;
        r[2] = rr[idxMicrophone + 2*nMicrophones] / cTs;

        n1 = (int) ceil(nSamples/(2*L[0]));
        n2 = (int) ceil(nSamples/(2*L[1]));
        n3 = (int) ceil(nSamples/(2*L[2]));

        // Generate room impulse response
        // Loop all possible duplicagtions on the x,y, and z axes:
        for (mx = -n1 ; mx <= n1 ; mx++)// x axis
        {
            Rm[0] = 2*mx*L[0];

            for (my = -n2 ; my <= n2 ; my++)// y axis
            {
                Rm[1] = 2*my*L[1];

                for (mz = -n3 ; mz <= n3 ; mz++)// z axis
                {
                    Rm[2] = 2*mz*L[2];
                    
                    // Loop all possible imaging variation for the axes:
                    for (q = 0 ; q <= 1 ; q++)// x axis
                    {
                        Rp_plus_Rm[0] = (1-2*q)*s[0] - r[0] + Rm[0];// x axis distance between virtual source and fixed reciver
                        xp[0] = 2*mx*LL_s[0] + (1-2*q)*ss[idxMicrophone + 0*nMicrophones];// virtual source x position
                        refl[0] = pow(beta[0], abs(mx-q)) * pow(beta[1], abs(mx));

                        for (j = 0 ; j <= 1 ; j++)// y axis
                        {
                            Rp_plus_Rm[1] = (1-2*j)*r[1] - s[1] + Rm[1];// y axis distance between virtual source and fixed reciver
                            xp[1] = 2*my*LL_s[1] + (1-2*j)*ss[idxMicrophone + 1*nMicrophones];// virtual source y position
                            refl[1] = pow(beta[2], abs(my-j)) * pow(beta[3], abs(my));

                            for (k = 0 ; k <= 1 ; k++)// z axis
                            {
                                Rp_plus_Rm[2] = (1-2*k)*r[2] - s[2] + Rm[2];// z axis distance between virtual source and fixed reciver
                                xp[2] = 2*mz*LL_s[2] + (1-2*k)*ss[idxMicrophone + 2*nMicrophones];// virtual source z position
                                refl[2] = pow(beta[4],abs(mz-k)) * pow(beta[5], abs(mz));
                                
                                order = abs(2*mx-q)+abs(2*my-j)+abs(2*mz-k);

                                dist = sqrt(pow(Rp_plus_Rm[0], 2) + pow(Rp_plus_Rm[1], 2) + pow(Rp_plus_Rm[2], 2));// Audio traveling distance
                                
                                face_num = box_ray(LL_s, xp, rr); // Determine the penetration face number
                                
                                if (face_num == 0)
                                    continue;// Do not proceed for artifact pathes
                                
                                // Translate coordinate system to receiver's room
                                r_rIM = r;
                                s_rIM = xp;
                                r_rIM[0] = r_rIM[0] - LL_s[0];
                                s_rIM[0] = s_rIM[0] - LL_s[0];
                                // Call receiver image method with the current source as a fixed source
                                ReciverIM(c, fs, rr, nMicrophones, ss, LL_r, beta_input_r, nDimension, nOrder, microphone_type, nSamples, angle, imp);
                                
                                // Calculate the sources's room reflections
                                if (abs(2*mx-q)+abs(2*my-j)+abs(2*mz-k) <= nOrder || nOrder == -1)
                                {
                                    fdist = floor(dist);
                                    
                                    if (fdist < nSamples)
                                    {
                                        // Calculate gain with respect to traveling distance and reflection order
                                        gain = sim_microphone(Rp_plus_Rm[0], Rp_plus_Rm[1], Rp_plus_Rm[2], angle, microphone_type[0])
                                            * refl[0]*refl[1]*refl[2]/(4*M_PI*dist*cTs);
                                        
                                        for (n = 0 ; n < Tw ; n++)
                                            LPI[n] =  0.5 * (1 - cos(2*M_PI*((n+1-(dist-fdist))/Tw))) * Fc * sinc(M_PI*Fc*(n+1-(dist-fdist)-(Tw/2)));
                                        
                                        startPosition = (int) fdist-(Tw/2)+1;
                                        for (n = 0 ; n < Tw; n++){
                                            if (startPosition+n >= 0 && startPosition+n < nSamples){
                                                imp[idxMicrophone + nMicrophones*(startPosition+n)] += gain * LPI[n];
                                            }
                                            
                                        }// n for loop
                                    }// conditioning for number of samples
                                }// conditioning for none-artifact sources
                            }// k for loop
                        }// j for loop
                    }// q for loop
                }// mz for loop
            }// my for loop
        }// mx for loop
        
        // 'Original' high-pass filter as proposed by Allen and Berkley.
        if (isHighPassFilter == 1)
        {
            for (int idx = 0 ; idx < 3 ; idx++) {Y[idx] = 0;}            
            for (int idx = 0 ; idx < nSamples ; idx++)
            {
                X0 = imp[idxMicrophone+nMicrophones*idx];
                Y[2] = Y[1];
                Y[1] = Y[0];
                Y[0] = B1*Y[1] + B2*Y[2] + X0;
                imp[idxMicrophone+nMicrophones*idx] = Y[0] + A1*Y[1] + R1*Y[2];
            }
        }
    }// Microphone for loop
    
    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        double* beta_hat = mxGetPr(plhs[1]);
        if (reverberation_time != 0) {
            beta_hat[0] = beta[0];
        }
        else {
            beta_hat[0] = 0;
        }
    }
    
    delete[] beta;
    delete[] microphone_type;
    delete[] Y;
    delete[] LPI;            
    delete[] r;
    delete[] s;
    delete[] L;    
}

