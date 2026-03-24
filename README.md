# EXP 3 : IIR-CHEBYSHEV-FITER-DESIGN

## AIM: 

 To design an IIR Chebyshev filter  using SCILAB. 

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM (LPF): 
```
clc;
clear;
close;

// ================== USER INPUT ==================
wp     = input('Enter the pass band frequency (Radians)= ');
ws     = input('Enter the stop band frequency (Radians)= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T      = input('Enter the value of sampling time= ');

// ================== PRE-WARPING ==================
omegap = (2/T)*tan(wp/2);
omegas = (2/T)*tan(ws/2);

// ================== FILTER ORDER ==================
Epsilon = sqrt((10^(0.1*alphap))-1);

N_calc = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) ...
         / acosh(omegas/omegap);

N = ceil(N_calc);

// ================== CUTOFF FREQUENCY ==================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));

// ================== ANALOG LOW PASS PROTOTYPE ==================
[pols, gn] = zpch1(N, Epsilon, omegac);
hs = syslin('c', gn, real(poly(pols, 's')));

// ================== BILINEAR TRANSFORMATION ==================
z = poly(0, 'z');
Hz = horner(hs, (2/T)*((1 - z^-1)/(1 + z^-1)));
Hz_sys = syslin('d', Hz);

// ================== NEAT CONSOLE OUTPUT ==================
mprintf("\n====================================================\n");
mprintf("      CHEBYSHEV TYPE-I DIGITAL LOW PASS FILTER\n");
mprintf("====================================================\n");

mprintf("\nInput Specifications:\n");
mprintf("----------------------------------------------------\n");
mprintf("Passband Frequency (wp)      = %8.4f rad/sample\n", wp);
mprintf("Stopband Frequency (ws)      = %8.4f rad/sample\n", ws);
mprintf("Passband Attenuation (Ap)    = %8.2f dB\n", alphap);
mprintf("Stopband Attenuation (As)    = %8.2f dB\n", alphas);
mprintf("Sampling Time (T)            = %8.6f sec\n", T);

mprintf("\nPrewarped Analog Frequencies:\n");
mprintf("----------------------------------------------------\n");
mprintf("Omega_p                      = %10.4f rad/sec\n", omegap);
mprintf("Omega_s                      = %10.4f rad/sec\n", omegas);

mprintf("\nFilter Order:\n");
mprintf("----------------------------------------------------\n");
mprintf("Calculated Order (N)         = %8.4f\n", N_calc);
mprintf("Rounded Filter Order         = %d\n", N);

mprintf("\nDesign Parameters:\n");
mprintf("----------------------------------------------------\n");
mprintf("Ripple Factor (Epsilon)      = %8.6f\n", Epsilon);
mprintf("Cutoff Frequency (Omega_c)   = %10.4f rad/sec\n", omegac);
mprintf("Prototype Gain               = %10.4e\n", gn);

disp("Digital Transfer Function H(z):");
disp(Hz_sys);

// ================== FREQUENCY RESPONSE ==================
[HW, w] = frmag(Hz_sys, 512);

figure();
plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency (w/pi)');
ylabel('Magnitude');
title('Frequency Response of Chebyshev Type-I Low Pass Filter');
```

## PROGRAM (HPF): 
```
clc;
clear;
close;

// ================== USER INPUT ==================
wp     = input('Enter the pass band frequency (Radians)= ');
ws     = input('Enter the stop band frequency (Radians)= ');
alphap = input('Enter the pass band attenuation (dB)= ');
alphas = input('Enter the stop band attenuation (dB)= ');
T      = input('Enter the value of sampling time= ');

// ================== PRE-WARPING ==================
omegap = (2/T)*tan(wp/2);
omegas = (2/T)*tan(ws/2);

// ================== FILTER ORDER ==================
Epsilon = sqrt((10^(0.1*alphap))-1);

N_calc = acosh(sqrt(((10^(0.1*alphas))-1)/((10^(0.1*alphap))-1))) ...
         / acosh(omegap/omegas);

N = ceil(N_calc);

// ================== CUTOFF FREQUENCY ==================
omegac = omegap / (((10^(0.1*alphap))-1)^(1/(2*N)));

// ================== ANALOG LOW PASS PROTOTYPE ==================
[pols, gn] = zpch1(N, Epsilon, omegac);
hs = syslin('c', gn, real(poly(pols, 's')));

// ================== LPF → HPF TRANSFORMATION ==================
s = poly(0, 's');
hs_hpf = horner(hs, (omegac^2 / s));

// ================== BILINEAR TRANSFORMATION ==================
z = poly(0, 'z');
Hz = horner(hs_hpf, (2/T)*((1 - z^-1)/(1 + z^-1)));
Hz_sys = syslin('d', Hz);

// ================== NEAT CONSOLE OUTPUT ==================
mprintf("\n====================================================\n");
mprintf("      CHEBYSHEV TYPE-I DIGITAL HIGH PASS FILTER\n");
mprintf("====================================================\n");

mprintf("\nInput Specifications:\n");
mprintf("----------------------------------------------------\n");
mprintf("Passband Frequency (wp)      = %8.4f rad/sample\n", wp);
mprintf("Stopband Frequency (ws)      = %8.4f rad/sample\n", ws);
mprintf("Passband Attenuation (Ap)    = %8.2f dB\n", alphap);
mprintf("Stopband Attenuation (As)    = %8.2f dB\n", alphas);
mprintf("Sampling Time (T)            = %8.6f sec\n", T);

mprintf("\nPrewarped Analog Frequencies:\n");
mprintf("----------------------------------------------------\n");
mprintf("Omega_p                      = %10.4f rad/sec\n", omegap);
mprintf("Omega_s                      = %10.4f rad/sec\n", omegas);

mprintf("\nFilter Order:\n");
mprintf("----------------------------------------------------\n");
mprintf("Calculated Order (N)         = %8.4f\n", N_calc);
mprintf("Rounded Filter Order         = %d\n", N);

mprintf("\nDesign Parameters:\n");
mprintf("----------------------------------------------------\n");
mprintf("Ripple Factor (Epsilon)      = %8.6f\n", Epsilon);
mprintf("Cutoff Frequency (Omega_c)   = %10.4f rad/sec\n", omegac);
mprintf("Prototype Gain               = %10.4e\n", gn);

// Optional display of transfer functions
disp("Analog High Pass Transfer Function:");
disp(hs_hpf);

disp("Digital Transfer Function H(z):");
disp(Hz_sys);

// ================== FREQUENCY RESPONSE ==================
[HW, w] = frmag(Hz_sys, 512);

figure();
plot(w/%pi, abs(HW));
xlabel('Normalized Digital Frequency (w/pi)');
ylabel('Magnitude');
title('Frequency Response of Chebyshev Type-I High Pass Filter');
```


## OUTPUT (LPF) : 
<img width="1915" height="1199" alt="image" src="https://github.com/user-attachments/assets/8fdfed81-098e-4251-bb18-888127690455" />

## OUTPUT (HPF) : 
<img width="1919" height="1080" alt="image" src="https://github.com/user-attachments/assets/368f3201-1d86-406a-b0d5-2499a3f318c7" />

## RESULT: 
The design of an IIR Chebyshev filter using SCILAB is sucessfully completed.
