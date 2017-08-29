close all
selectBackgroundCell = 1; %% Set = 1 for images where you want to subtract background and Set = 0 for images without backgroung

MHRun_2G_I4.Location = 'E:/Images/MH+DZ/03-09-2017/MHRun_2G_I4.czi'; 
MHRun_2G_I4 = Run(MHRun_2G_I4, selectBackgroundCell);

MHRun_11G_I4.Location = 'E:/Images/MH+DZ/03-09-2017/MHRun_11G_I4.czi'; 
MHRun_11G_I4 = Run(MHRun_11G_I4, selectBackgroundCell);

MHRun_11G_3MH_I4.Location = 'E:/Images/MH+DZ/03-09-2017/MHRun_11G_3MH_I4.czi'; 
MHRun_11G_3MH_I4 = Run(MHRun_11G_3MH_I4, selectBackgroundCell);

%% Activity Calculations

DZRun_data=[
    DZRun_2G_I1.RatioActive DZRun_11G_I1.RatioActive DZRun_11G_50DZ_I1.RatioActive;
    DZRun_2G_I2.RatioActive DZRun_11G_I2.RatioActive DZRun_11G_50DZ_I2.RatioActive;
    DZRun_2G_I3.RatioActive DZRun_11G_I3.RatioActive DZRun_11G_50DZ_I3.RatioActive; 
    DZRun_2G_I4.RatioActive DZRun_11G_I4.RatioActive DZRun_11G_50DZ_I4.RatioActive; 
]

MHRun_data=[
    MHRun_2G_I1.RatioActive MHRun_11G_I1.RatioActive MHRun_11G_3MH_I1.RatioActive;
    MHRun_2G_I2.RatioActive MHRun_11G_I2.RatioActive MHRun_11G_3MH_I2.RatioActive;
    MHRun_2G_I3.RatioActive MHRun_11G_I3.RatioActive MHRun_11G_3MH_I3.RatioActive; 
    MHRun_2G_I4.RatioActive MHRun_11G_I4.RatioActive MHRun_11G_3MH_I4.RatioActive; 
]

DZMHRun_data=[
    DZMHRun_2G_I1.RatioActive DZMHRun_11G_I1.RatioActive DZMHRun_11G_50DZ_3MH_I1.RatioActive;
    DZMHRun_2G_I2.RatioActive DZMHRun_11G_I2.RatioActive DZMHRun_11G_50DZ_3MH_I2.RatioActive;
    DZMHRun_2G_I3.RatioActive DZMHRun_11G_I3.RatioActive DZMHRun_11G_50DZ_3MH_I3.RatioActive; 
    DZMHRun_2G_I5.RatioActive DZMHRun_11G_I5.RatioActive DZMHRun_11G_50DZ_3MH_I5.RatioActive; 
    DZMHRun_2G_I6.RatioActive DZMHRun_11G_I6.RatioActive DZMHRun_11G_50DZ_3MH_I6.RatioActive; 
]


%% Correlation Calculations 

DZRun_data=[
    DZRun_2G_I1.RatioCorr_3 DZRun_11G_I1.RatioCorr_3 DZRun_11G_50DZ_I1.RatioCorr_3;
    DZRun_2G_I2.RatioCorr_3 DZRun_11G_I2.RatioCorr_3 DZRun_11G_50DZ_I2.RatioCorr_3;
    DZRun_2G_I3.RatioCorr_3 DZRun_11G_I3.RatioCorr_3 DZRun_11G_50DZ_I3.RatioCorr_3; 
    DZRun_2G_I4.RatioCorr_3 DZRun_11G_I4.RatioCorr_3 DZRun_11G_50DZ_I4.RatioCorr_3; 
]

MHRun_data=[
    MHRun_2G_I1.RatioCorr_3 MHRun_11G_I1.RatioCorr_3 MHRun_11G_3MH_I1.RatioCorr_3;
    MHRun_2G_I2.RatioCorr_3 MHRun_11G_I2.RatioCorr_3 MHRun_11G_3MH_I2.RatioCorr_3;
    MHRun_2G_I3.RatioCorr_3 MHRun_11G_I3.RatioCorr_3 MHRun_11G_3MH_I3.RatioCorr_3; 
    MHRun_2G_I4.RatioCorr_3 MHRun_11G_I4.RatioCorr_3 MHRun_11G_3MH_I4.RatioCorr_3; 
]

DZMHRun_data=[
    DZMHRun_2G_I1.RatioCorr_3 DZMHRun_11G_I1.RatioCorr_3 DZMHRun_11G_50DZ_3MH_I1.RatioCorr_3;
    DZMHRun_2G_I2.RatioCorr_3 DZMHRun_11G_I2.RatioCorr_3 DZMHRun_11G_50DZ_3MH_I2.RatioCorr_3;
    DZMHRun_2G_I3.RatioCorr_3 DZMHRun_11G_I3.RatioCorr_3 DZMHRun_11G_50DZ_3MH_I3.RatioCorr_3; 
    DZMHRun_2G_I5.RatioCorr_3 DZMHRun_11G_I5.RatioCorr_3 DZMHRun_11G_50DZ_3MH_I5.RatioCorr_3; 
    DZMHRun_2G_I6.RatioCorr_3 DZMHRun_11G_I6.RatioCorr_3 DZMHRun_11G_50DZ_3MH_I6.RatioCorr_3; 
]



