.. highlight:: matlab
.. _example_generatesignals:

*****************************************
Generating dipolar signals
*****************************************


.. code-block:: matlab

    %======================================================================
    % DeerAnalyis2
    % Example: Generating dipolar signals
    % How to generate dipolar signals corresponding to different experiments
    %=======================================================================

    clear,clc,clf

    %Preparation
    %----------------------------------------------
    N = 200;
    t = linspace(-0.2,8,N);
    r = time2dist(t);
    P = rd_onegaussian(r,[4,0.3]);

    %4-pulse DEER
    %----------------------------------------------
    V4pdeer = dipolarsignal(t,r,P,'ModDepth',0.5);

    %5-pulse DEER
    %----------------------------------------------
    %Simulate different artefact levels
    V5pdeer1 = dipolarsignal(t,r,P,'ModDepth',0.5,'FivePulseCoeff',0.3);
    V5pdeer2 = dipolarsignal(t,r,P,'ModDepth',0.5,'FivePulseCoeff',0.2);
    V5pdeer3 = dipolarsignal(t,r,P,'ModDepth',0.5,'FivePulseCoeff',0.1);

    %RIDME
    %----------------------------------------------
    Tmix = 50; %Mixing time [us]
    T1 = 88; %Relaxation time [us]
    %Simulate different overtones
    OverCoeff = overtones(1,Tmix,T1);
    Vridme1 = dipolarsignal(t,r,P,'ModDepth',0.5,'Overtones',OverCoeff);
    OverCoeff = overtones(2,Tmix,T1);
    Vridme2 = dipolarsignal(t,r,P,'ModDepth',0.5,'Overtones',OverCoeff);
    OverCoeff = overtones(3,Tmix,T1);
    Vridme3 = dipolarsignal(t,r,P,'ModDepth',0.5,'Overtones',OverCoeff);


    %Plotting
    %----------------------------------------------

    subplot(141)
    plot(r,P,'k','LineWidth',1)
    axis tight, box on, grid on
    title('DIstance Distribution')
    xlabel('Distance [nm]')
    ylabel('P(r)')

    subplot(142)
    plot(t,V4pdeer,'LineWidth',1)
    axis tight, box on, grid on
    title('4-pulse DEER')
    xlabel('Time [\mus]')
    ylabel('V(t)')

    subplot(143)
    plot(t,V5pdeer1,t,V5pdeer2,t,V5pdeer3,'LineWidth',1)
    axis tight, box on, grid on
    legend('\rho = 0.3','\rho = 0.2','\rho = 0.1')
    title('5-pulse DEER')
    xlabel('Time [\mus]')
    ylabel('V(t)')

    subplot(144)
    plot(t,Vridme1,t,Vridme2,t,Vridme3,'LineWidth',1)
    axis tight, box on, grid on
    legend('1 Overtone','2 Overtones','3 Overtones')
    title('RIDME')
    xlabel('Time [\mus]')
    ylabel('V(t)')

.. figure:: ../images/example_generate_signals.svg
    :align: center