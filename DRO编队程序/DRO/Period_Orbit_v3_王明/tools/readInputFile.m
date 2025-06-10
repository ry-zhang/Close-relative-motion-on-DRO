
function aux = readInputFile(aux)
%
% ���ܣ���excel�ļ���ȡ������Ͳ���
%
% ���룺inputFile.xlsx�ļ�
%
% �����aux�ṹ��
%
% ���ߣ�����
% ��λ���п�Ժ�ռ�Ӧ�ù����뼼������
% ʱ�䣺2021��1��17��
%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------ �����������ļ� -------------------
addpath('tools');
addpath('eph');
load('DE430.mat')
aux.DE430 = DE430;
try
    sheet_solarSys = readtable('inputFile.xlsx' , 'Sheet' , 'solarSys' , 'Format' , 'auto' , 'ReadRowNames' , true);
    sheet_aux = readtable('inputFile.xlsx' , 'Sheet' , 'aux' , 'Format' , 'auto' , 'ReadRowNames' , true);
    sheet_PeriOrbit = readtable('inputFile.xlsx' , 'Sheet' , 'PeriOrbit' , 'Format' , 'auto' , 'ReadRowNames' , true);
    sheet_PeriOribit_Modi = readtable('inputFile.xlsx' , 'Sheet' , 'PeriOrbit_Modi' , 'Format' , 'auto' , 'ReadRowNames' , true);
catch
    % �������ʹ���������
    sheet_solarSys = readtable('inputFile.xlsx' , 'Sheet' , 'solarSys' , 'ReadRowNames' , true);
    sheet_aux = readtable('inputFile.xlsx' , 'Sheet' , 'aux' , 'ReadRowNames' , true);
    sheet_PeriOrbit = readtable('inputFile.xlsx' , 'Sheet' , 'PeriOrbit' , 'ReadRowNames' , true);
    sheet_PeriOribit_Modi = readtable('inputFile.xlsx' , 'Sheet' , 'PeriOrbit_Modi' , 'ReadRowNames' , true);
end



% --------------------- ̫��ϵ���� ----------------------
aux.planet.name = eval(sheet_solarSys({'name'},:).value{1});
aux.planet.mu = eval(sheet_solarSys({'mu'},:).value{1});
aux.planet.req = eval(sheet_solarSys({'req'},:).value{1});
aux.planet.soi = eval(sheet_solarSys({'soi'},:).value{1});
aux.planet.AU = eval(sheet_solarSys({'AU'},:).value{1});
aux.planet.g0 = eval(sheet_solarSys({'g0'},:).value{1});


% -------------------- ����ϵͳ���� --------------------
aux.LU = eval(sheet_aux({'LU'},:).value{1}); % (km)
aux.TU = eval(sheet_aux({'TU'},:).value{1}); % (s)
aux.VU = eval(sheet_aux({'VU'},:).value{1}); % (km/s)
aux.mu = eval(sheet_aux({'mu'},:).value{1}); % (--)
aux.ms = eval(sheet_aux({'ms'},:).value{1}); % (--)
aux.rhos = eval(sheet_aux({'rhos'},:).value{1}); % (--)
aux.oms = eval(sheet_aux({'oms'},:).value{1}); % (--)
aux.tol = eval(sheet_aux({'tol'},:).value{1}); % (--)



% ----------------------- ������� -----------------------
% �������
aux.orb_type = sheet_aux({'orb_type'},:).value{1};

%  ����ֲ�
aux.Burf_3D = eval(sheet_aux({'Burf_3D'},:).value{1});

% ���ش�нڵ���
aux.node_n = eval(sheet_aux({'node_n'},:).value{1});

% α��������
aux.delta_s = eval(sheet_aux({'delta_s'},:).value{1});

% α�������ش���
aux.conti_n = eval(sheet_aux({'conti_n'},:).value{1});

if aux.node_n < 2
    fprintf('����! aux.node_n �������2! \n');
    aux.err_IO = 1;
    prompt = '����������n��2��ֵ��������n:\n';
    n = input(prompt);
    if n>=2
        aux.node_n = n;
        aux.err_IO = 0;
    else
        fprintf('û��Ҫ�����룬���ź���������ֹ! \n');
       return;
      
    end
else
    aux.err_IO = 0;
end

% -------------- ���ڹ�� --------------
% ����DRO��ֵ (LU , VU)
aux.DRO_xx = eval(sheet_PeriOrbit({'DRO_xx'},:).value{1});

%����DRO���� (TU)
aux.DRO_P = eval(sheet_PeriOrbit({'DRO_P'},:).value{1});

%--------------DPO-------------------------
% DPO��ֵ (LU , VU)
aux.DPO_xx = eval(sheet_PeriOrbit({'DPO_xx'},:).value{1});

% DPO���� (TU)
aux.DPO_P = eval(sheet_PeriOrbit({'DPO_P'},:).value{1});

%-----------------LoPO---------------------
% LoPO��ֵ (LU , VU)
aux.LoPO_xx = eval(sheet_PeriOrbit({'LoPO_xx'},:).value{1});

% LoPO���� (TU)
aux.LoPO_P = eval(sheet_PeriOrbit({'LoPO_P'},:).value{1});

% -------------- Halo/Lyapunov��� --------------
% (/) Halo/Lyapunov��� ƽ��������

aux.Halo_type = sheet_aux({'Halo_type'},:).value{1};
aux.Halo_class = sheet_aux({'Halo_class'},:).value{1};
aux.Lyapunov_type = sheet_aux({'Lyapunov_type'},:).value{1};
aux.Lyapunov_xx = eval(sheet_PeriOrbit({'Lyapunov_xx'},:).value{1});
aux.Lyapunov_P = eval(sheet_PeriOrbit({'Lyapunov_P'},:).value{1});


% ---------------------������ ----------------------------
% ����Ȧ��(/)
aux.resoOrb_p = eval(sheet_aux({'resoOrb_p'},:).value{1});
% С����Ȧ��(/)
aux.resoOrb_q = eval(sheet_aux({'resoOrb_q'},:).value{1});

%--------------------- ���������г�ʼ���� ----------------------
temp_resoOrb_data = eval(sheet_PeriOrbit({'resoOrb_type'},:).value{1});

%---------------------���������Ϳ�---------------------------
temp_resoOrb_type = temp_resoOrb_data (:,1);

%---------------------���������ص�߶ȿ�(km)---------------------------
temp_resoOrb_rp = temp_resoOrb_data (:,2);

%---------------------��p��q�����Լ��---------------------------
[~, gcd] = is_coprime(aux.resoOrb_p  , aux.resoOrb_q);

if aux.resoOrb_p/gcd>4 || aux.resoOrb_q>4
    fprintf('--------------------------------------------------- \n');
    fprintf('ע��! Ŀǰ��֧��p��q��Լ���С��4����� \n');
    aux.err_IO = 1;
    prompt = '��������������p��������p:\n';
    p = input(prompt);
    prompt = '��������������q��������q:\n';
    q = input(prompt);
    [~, gcd] = is_coprime(p  , q);
    if  p/gcd<=4 && q/gcd<=4
        aux.resoOrb_p = p;
        aux.resoOrb_q = q;
        aux.err_IO = 0;
    else
         fprintf('û��Ҫ�����룬���ź���������ֹ! \n');
       return;       
    end     
end
% ����������
aux.resoOrb_type = "M"+num2str(aux.resoOrb_p/gcd)+"N"+num2str(aux.resoOrb_q/gcd);
aux.gcd = gcd;
% ƥ��������
index = strcmp(aux.resoOrb_type,temp_resoOrb_type);

% ���ص�߶�(LU)
aux.resoOrb_rp = str2double(temp_resoOrb_rp(index))/aux.LU;

% ������(rad)
aux.resoOrb_incl = eval(sheet_aux({'resoOrb_incl'},:).value{1}) * pi / 180;

% ������(rad)
% aux.resoOrb_raan = eval(sheet_aux({'resoOrb_raan'},:).value{1}) * pi / 180;
aux.resoOrb_raan = 0;

% ������(rad)
aux.resoOrb_tanom = eval(sheet_aux({'resoOrb_tanom'},:).value{1}) * pi / 180;

% ---------------------------------- ����ƽ���� ---------------------------------
aux.TriangleL_type = sheet_aux({'TriangleL_type'},:).value{1};
aux.TriangleL_class = sheet_aux({'TriangleL_class'},:).value{1};
aux.TriangleL_xx = eval(sheet_PeriOrbit({'TriangleL_xx'},:).value{1});
aux.TriangleL_P = eval(sheet_PeriOrbit({'TriangleL_P'},:).value{1});

% ---------------------- NLP���� ------------------------
aux.nlp.algorithm = sheet_aux({'nlp_algorithm'},:).value{1};
aux.nlp.maxIter = eval(sheet_aux({'nlp_maxIter'},:).value{1});
aux.nlp.maxFunEvals = eval(sheet_aux({'nlp_maxFunEvals'},:).value{1});
aux.nlp.tolX = eval(sheet_aux({'nlp_tolX'},:).value{1});
aux.nlp.tolCon = eval(sheet_aux({'nlp_tolCon'},:).value{1});
aux.nlp.tolFun = eval(sheet_aux({'nlp_tolFun'},:).value{1});
aux.method = sheet_aux({'method'},:).value{1};
aux.grad_IO = sheet_aux({'grad_IO'},:).value{1};


% ---------------------- ���ڹ������/��������------------------------
% �������������
aux.periOrb_Modi =  sheet_PeriOribit_Modi({'periOrb_Modi'},:).value{1};

aux.periOrb_T = sheet_PeriOribit_Modi({'periOrb_T'},:).value{1};
% p,q�����
aux.periOrb_p = eval(sheet_PeriOribit_Modi({'periOrb_p'},:).value{1});
aux.periOrb_q = eval(sheet_PeriOribit_Modi({'periOrb_q'},:).value{1});
% m, �������ڱ���
aux.periOrb_m = eval(sheet_PeriOribit_Modi({'periOrb_m'},:).value{1});
% n, ÿȦ����ڵ���
aux.periOrb_n = eval(sheet_PeriOribit_Modi({'periOrb_n'},:).value{1});
aux.Sun_angle = eval(sheet_PeriOribit_Modi({'Sun_angle'},:).value{1});
aux.Sun_angle = deg2rad(aux.Sun_angle);                                  % �Ƕ�ת����                                     
aux.dep_t_jdtdb = eval(sheet_PeriOribit_Modi({'dep_t_jdtdb'},:).value{1});   % ��׼ʱ��
end
