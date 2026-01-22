% 该代码为生成各种扫描类型的电场强度 为后续计算雅可比矩阵
% function out = model
clear all
clc
tic
N = 1;
for n = 1:N
%%
% Define the conductivities of the objects and the conductivity of
% background
% cond1 = [1];
% cond2 = [0.1];
cond_background = [2];

disp('---------------------------------------------');

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('F:\Desk_left\EIT_MATLAB\1_Array\ref');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

%创建几何
model.component('comp1').geom('geom1').useConstrDim(false);
model.component('comp1').geom('geom1').repairTolType('relative');
model.component('comp1').geom('geom1').create('imp1', 'Import');
model.component('comp1').geom('geom1').feature('imp1').set('filename', 'F:\Desk_left\EIT_MATLAB\1_Array\EMP16.mphbin');
model.component('comp1').geom('geom1').feature('imp1').importData;
model.component('comp1').geom('geom1').run('imp1');

%创建材料1-电极材料参数
model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').label('Titanium beta-21S');
model.component('comp1').material('mat1').set('family', 'titanium');
model.component('comp1').material('mat1').propertyGroup('def').label('Basic');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').descr('relpermeability_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'7.407e5[S/m]' '0' '0' '0' '7.407e5[S/m]' '0' '0' '0' '7.407e5[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').descr('electricconductivity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'7.06e-6[1/K]' '0' '0' '0' '7.06e-6[1/K]' '0' '0' '0' '7.06e-6[1/K]'});
model.component('comp1').material('mat1').propertyGroup('def').descr('thermalexpansioncoefficient_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '710[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').descr('heatcapacity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').descr('relpermittivity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '4940[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('def').descr('density_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'7.5[W/(m*K)]' '0' '0' '0' '7.5[W/(m*K)]' '0' '0' '0' '7.5[W/(m*K)]'});
model.component('comp1').material('mat1').propertyGroup('def').descr('thermalconductivity_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('Enu').label('Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').propertyGroup('Enu').set('youngsmodulus', '105e9[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').descr('youngsmodulus_symmetry', '');
model.component('comp1').material('mat1').propertyGroup('Enu').set('poissonsratio', '0.33');
model.component('comp1').material('mat1').propertyGroup('Enu').descr('poissonsratio_symmetry', '');
model.component('comp1').material('mat1').set('groups', {});
model.component('comp1').material('mat1').set('family', 'titanium');
model.component('comp1').geom('geom1').run;
model.component('comp1').material('mat1').selection.set([1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);

%创建材料2—参考水
model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat2').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat2').propertyGroup('def').func.create('rho', 'Piecewise');
model.component('comp1').material('mat2').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat2').propertyGroup('def').func.create('cs', 'Interpolation');
model.component('comp1').material('mat2').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat2').propertyGroup('def').func.create('an2', 'Analytic');
model.component('comp1').material('mat2').propertyGroup('def').func.create('an3', 'Analytic');
model.component('comp1').material('mat2').label('Water, liquid');
model.component('comp1').material('mat2').set('family', 'water');
model.component('comp1').material('mat2').propertyGroup('def').label('Basic');
model.component('comp1').material('mat2').propertyGroup('def').func('eta').label('Piecewise');
model.component('comp1').material('mat2').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat2').propertyGroup('def').func('eta').set('pieces', {'273.15' '413.15' '1.3799566804-0.021224019151*T^1+1.3604562827E-4*T^2-4.6454090319E-7*T^3+8.9042735735E-10*T^4-9.0790692686E-13*T^5+3.8457331488E-16*T^6'; '413.15' '553.75' '0.00401235783-2.10746715E-5*T^1+3.85772275E-8*T^2-2.39730284E-11*T^3'});
model.component('comp1').material('mat2').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat2').propertyGroup('def').func('Cp').label('Piecewise 2');
model.component('comp1').material('mat2').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat2').propertyGroup('def').func('Cp').set('pieces', {'273.15' '553.75' '12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4'});
model.component('comp1').material('mat2').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat2').propertyGroup('def').func('rho').label('Piecewise 3');
model.component('comp1').material('mat2').propertyGroup('def').func('rho').set('arg', 'T');
model.component('comp1').material('mat2').propertyGroup('def').func('rho').set('smooth', 'contd1');
model.component('comp1').material('mat2').propertyGroup('def').func('rho').set('pieces', {'273.15' '293.15' '0.000063092789034*T^3-0.060367639882855*T^2+18.9229382407066*T-950.704055329848'; '293.15' '373.15' '0.000010335053319*T^3-0.013395065634452*T^2+4.969288832655160*T+432.257114008512'});
model.component('comp1').material('mat2').propertyGroup('def').func('rho').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat2').propertyGroup('def').func('k').label('Piecewise 4');
model.component('comp1').material('mat2').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat2').propertyGroup('def').func('k').set('pieces', {'273.15' '553.75' '-0.869083936+0.00894880345*T^1-1.58366345E-5*T^2+7.97543259E-9*T^3'});
model.component('comp1').material('mat2').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat2').propertyGroup('def').func('cs').label('Interpolation');
model.component('comp1').material('mat2').propertyGroup('def').func('cs').set('table', {'273' '1403';  ...
'278' '1427';  ...
'283' '1447';  ...
'293' '1481';  ...
'303' '1507';  ...
'313' '1526';  ...
'323' '1541';  ...
'333' '1552';  ...
'343' '1555';  ...
'353' '1555';  ...
'363' '1550';  ...
'373' '1543'});
model.component('comp1').material('mat2').propertyGroup('def').func('cs').set('interp', 'piecewisecubic');
model.component('comp1').material('mat2').propertyGroup('def').func('cs').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat2').propertyGroup('def').func('an1').label('Analytic 1');
model.component('comp1').material('mat2').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat2').propertyGroup('def').func('an1').set('expr', '-1/rho(T)*d(rho(T),T)');
model.component('comp1').material('mat2').propertyGroup('def').func('an1').set('args', {'T'});
model.component('comp1').material('mat2').propertyGroup('def').func('an1').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat2').propertyGroup('def').func('an1').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat2').propertyGroup('def').func('an2').label('Analytic 2');
model.component('comp1').material('mat2').propertyGroup('def').func('an2').set('funcname', 'gamma_w');
model.component('comp1').material('mat2').propertyGroup('def').func('an2').set('expr', '1+(T/Cp(T))*(alpha_p(T)*cs(T))^2');
model.component('comp1').material('mat2').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat2').propertyGroup('def').func('an2').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('an2').set('fununit', '1');
model.component('comp1').material('mat2').propertyGroup('def').func('an2').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat2').propertyGroup('def').func('an3').label('Analytic 3');
model.component('comp1').material('mat2').propertyGroup('def').func('an3').set('funcname', 'muB');
model.component('comp1').material('mat2').propertyGroup('def').func('an3').set('expr', '2.79*eta(T)');
model.component('comp1').material('mat2').propertyGroup('def').func('an3').set('args', {'T'});
model.component('comp1').material('mat2').propertyGroup('def').func('an3').set('argunit', 'K');
model.component('comp1').material('mat2').propertyGroup('def').func('an3').set('fununit', 'Pa*s');
model.component('comp1').material('mat2').propertyGroup('def').func('an3').set('plotargs', {'T' '273.15' '553.75'});
model.component('comp1').material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat2').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(T)' '0' '0' '0' 'alpha_p(T)' '0' '0' '0' 'alpha_p(T)'});
model.component('comp1').material('mat2').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat2').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat2').propertyGroup('def').descr('dynamicviscosity_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').set('ratioofspecificheat', 'gamma_w(T)');
model.component('comp1').material('mat2').propertyGroup('def').descr('ratioofspecificheat_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]' '0' '0' '0' '5.5e-6[S/m]'});
model.component('comp1').material('mat2').propertyGroup('def').descr('electricconductivity_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat2').propertyGroup('def').descr('heatcapacity_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').set('density', 'rho(T)');
model.component('comp1').material('mat2').propertyGroup('def').descr('density_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
model.component('comp1').material('mat2').propertyGroup('def').descr('thermalconductivity_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat2').propertyGroup('def').descr('soundspeed_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').descr('thermalexpansioncoefficient_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').descr('bulkviscosity_symmetry', '');
model.component('comp1').material('mat2').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat2').set('groups', {});
model.component('comp1').material('mat2').set('family', 'water');
model.component('comp1').material('mat2').materialType('nonSolid');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'80'});
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {strcat(num2str(cond_background(1)),'[S/m]')});
model.component('comp1').material('mat2').selection.set([17]);

%设置物理场-电流
model.component('comp1').physics.create('ec', 'ConductiveMedia', 'geom1');%创建

model.component('comp1').physics('ec').create('gnd1', 'Ground', 0);
model.component('comp1').physics('ec').feature('gnd1').selection.set([41]);%接地端的点 编号41-由comsol软件建立后自动排的

model.component('comp1').physics('ec').create('term1', 'Terminal', 1);
model.component('comp1').physics('ec').feature('term1').selection.set([25]);%正向输入电流-线段-编号25
model.component('comp1').physics('ec').feature('term1').set('I0', 1);%设置输入电流1

model.component('comp1').physics('ec').create('term2', 'Terminal', 1);
model.component('comp1').physics('ec').feature('term2').selection.set([31]);%负向输入电流-线段-编号31
model.component('comp1').physics('ec').feature('term2').set('I0', -1);%设置输入电流-1

%设置自动网格
% model.component('comp1').mesh('mesh1').autoMeshSize(1); %1是极细 2是
% model.component('comp1').mesh('mesh1').run;

%设置自定义网格
model.component('comp1').mesh('mesh1').automatic(false);
model.component('comp1').mesh('mesh1').feature('size').set('custom', false);
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', true);
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.0054);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 4.09E-6);
model.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1);
model.component('comp1').mesh('mesh1').feature('size').set('hcurve', 0.1);
model.component('comp1').mesh('mesh1').run;


%%
% Create solve first 
model.study.create('std1');
model.study('std1').create('stat', 'Stationary');%创建稳态
model.study('std1').feature('stat').activate('ec', true);

model.sol.create('sol1');%解
model.result.create('pg1', 'PlotGroup2D');%创建2D绘图
model.result('pg1').label([native2unicode(hex2dec({'57' '3a'}), 'unicode')  native2unicode(hex2dec({'5f' '3a'}), 'unicode') ]);
model.result('pg1').feature.create('surf1', 'Surface');%表面绘图

% Create the electrode index
eNum=16;
%该编号对应电极1-16-1的边的编号，之所以是下列编号是因为comsol绘图自动排序的
indexMap=[25 31 37 43 48 42 36 30 23 17 11 5 1 7 13 19 25]; % The last one is the first one to form a loop
cIndex=1;

for i=1:eNum
% for i=1:eNum
    % Select the excited electrode
    model.physics('ec').feature('term1').selection.set([]); % select current 1
    model.physics('ec').feature('term2').selection.set([]); % select current 2
    
    model.physics('ec').feature('term1').selection.set([indexMap(i)]);   % select current 1
    model.physics('ec').feature('term2').selection.set([indexMap(i+1)]); % select current 2

    model.sol('sol1').study('std1');
    model.sol('sol1').feature.remove('s1');  % Causion
    model.sol('sol1').feature.remove('v1');  % Causion
    model.sol('sol1').feature.remove('st1'); % Causion

    model.study('std1').feature('stat').set('notlistsolnum', 1);
    model.study('std1').feature('stat').set('notsolnum', '1');
    model.study('std1').feature('stat').set('listsolnum', 1);
    model.study('std1').feature('stat').set('solnum', '1');

    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').feature('st1').set('study', 'std1');
    model.sol('sol1').feature('st1').set('studystep', 'stat');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').feature('v1').set('control', 'stat');
    model.sol('sol1').create('s1', 'Stationary');
    model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
    model.sol('sol1').feature('s1').feature.remove('fcDef');
    model.sol('sol1').attach('std1');

    model.result('pg1').feature('surf1').set('data', 'dset1');
    model.result('pg1').feature('surf1').set('expr', 'ec.normE');%设置电场模
    model.result('pg1').feature('surf1').set('descr', [native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode')  native2unicode(hex2dec({'6a' '21'}), 'unicode') ]);


    model.sol('sol1').runAll;
    model.result('pg1').run;
    
    figure;
    mphplot(model,'pg1','rangenum',1);%函数在说明手册314页
%     axis square
    axis off
    
    if i>1
        T=eNum;
    else
        T=eNum-1;
    end
   
    % Measurement
    model.result.export.create('data1', 'Data');
    model.result.export('data1').set('expr', {'ec.normE'});
    model.result.export('data1').set('descr', {[native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode')  native2unicode(hex2dec({'6a' '21'}), 'unicode') ]});
    model.result.export('data1').set('unit', {'V/m'});
    model.result.export('data1').set('expr', {'ec.normE' 'ec.normE'});
    model.result.export('data1').set('descr', {[native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode')  native2unicode(hex2dec({'6a' '21'}), 'unicode') ] [native2unicode(hex2dec({'75' '35'}), 'unicode')  native2unicode(hex2dec({'57' '3a'}), 'unicode')  native2unicode(hex2dec({'6a' '21'}), 'unicode') ]});
    model.result.export('data1').remove('unit', 1);
    model.result.export('data1').remove('descr', 1);
    model.result.export('data1').remove('expr', [1]);
    model.result.export('data1').set('location', 'grid');
    model.result.export.remove('data1');
    model.result.export.create('tbl1', 'Table');
    model.result.export('tbl1').set('source', 'evaluationgroup');
    model.result.export.remove('tbl1');
    model.result('pg1').run;
    model.result.export.create('plot1', 'Plot');
    model.result.export('plot1').set('plotgroup', 'pg1');
    model.result.export('plot1').set('struct', 'spreadsheet');
    model.result.export('plot1').set('filename', strcat('.\data_E1\',sprintf('%04d',i),'.txt'));
%     model.result.export('plot1').set('filename', strcat('.\data_E\',num2str(i),'.txt'));
    model.result.export('plot1').run;
    model.result.export.remove('plot1');
    
end

%% Export the Measurement data



end
    