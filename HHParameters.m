classdef HHParameters
    properties
        C_m = 1.0     % 막 정전용량 (μF/cm²)
        
        % 전도도 (mS/cm²)
        g_Na = 120.0  % 나트륨
        g_K = 36.0    % 칼륨
        g_L = 0.3     % 누출
        
        % 평형전위 (mV)
        E_Na = 50.0   % 나트륨
        E_K = -77.0   % 칼륨
        E_L = -54.387 % 누출
    end
end