classdef HodgkinHuxleyModel < handle
    properties
        params
    end
    
    methods
        % 생성자
        function obj = HodgkinHuxleyModel()
            obj.params = HHParameters();
        end
        
        % 게이팅 변수의 속도 상수
        function val = alpha_m(obj, V)
            val = 0.1 * (V + 40) ./ (1 - exp(-(V + 40) / 10));
        end
        
        function val = beta_m(obj, V)
            val = 4.0 * exp(-(V + 65) / 18);
        end
        
        function val = alpha_h(obj, V)
            val = 0.07 * exp(-(V + 65) / 20);
        end
        
        function val = beta_h(obj, V)
            val = 1.0 ./ (1 + exp(-(V + 35) / 10));
        end
        
        function val = alpha_n(obj, V)
            val = 0.01 * (V + 55) ./ (1 - exp(-(V + 55) / 10));
        end
        
        function val = beta_n(obj, V)
            val = 0.125 * exp(-(V + 65) / 80);
        end
        
        % 정상상태 값 계산
        function [m_inf, h_inf, n_inf] = steady_state_values(obj, V)
            m_inf = obj.alpha_m(V) ./ (obj.alpha_m(V) + obj.beta_m(V));
            h_inf = obj.alpha_h(V) ./ (obj.alpha_h(V) + obj.beta_h(V));
            n_inf = obj.alpha_n(V) ./ (obj.alpha_n(V) + obj.beta_n(V));
        end
        
        % 이온 전류 계산
        function [I_Na, I_K, I_L] = ionic_currents(obj, V, m, h, n)
            I_Na = obj.params.g_Na * m.^3 .* h .* (V - obj.params.E_Na);
            I_K = obj.params.g_K * n.^4 .* (V - obj.params.E_K);
            I_L = obj.params.g_L * (V - obj.params.E_L);
        end
        
        % 미분방정식 시스템
        function [dVdt, dmdt, dhdt, dndt] = derivatives(obj, V, m, h, n, t, I_stim)
            [I_Na, I_K, I_L] = obj.ionic_currents(V, m, h, n);
            
            dVdt = (I_stim(t) - I_Na - I_K - I_L) / obj.params.C_m;
            dmdt = obj.alpha_m(V) .* (1 - m) - obj.beta_m(V) .* m;
            dhdt = obj.alpha_h(V) .* (1 - h) - obj.beta_h(V) .* h;
            dndt = obj.alpha_n(V) .* (1 - n) - obj.beta_n(V) .* n;
        end
        
        % 시뮬레이션 실행
        function [time, V, m, h, n] = simulate(obj, t_max, dt, V0, I_stim)
            if nargin < 2, t_max = 100; end
            if nargin < 3, dt = 0.01; end
            if nargin < 4, V0 = -65.0; end
            if nargin < 5
                I_stim = @(t) 10.0 * (t > 10 & t < 40);
            end
            
            time = 0:dt:t_max;
            N = length(time);
            V = zeros(1, N);
            m = zeros(1, N);
            h = zeros(1, N);
            n = zeros(1, N);
            
            V(1) = V0;
            [m(1), h(1), n(1)] = obj.steady_state_values(V0);
            
            for i = 1:N-1
                [dVdt, dmdt, dhdt, dndt] = obj.derivatives(V(i), m(i), h(i), n(i), time(i), I_stim);
                
                k1_V = dt * dVdt;
                k1_m = dt * dmdt;
                k1_h = dt * dhdt;
                k1_n = dt * dndt;
                
                [dVdt, dmdt, dhdt, dndt] = obj.derivatives(V(i) + k1_V/2, m(i) + k1_m/2, ...
                    h(i) + k1_h/2, n(i) + k1_n/2, time(i) + dt/2, I_stim);
                k2_V = dt * dVdt;
                k2_m = dt * dmdt;
                k2_h = dt * dhdt;
                k2_n = dt * dndt;
                
                [dVdt, dmdt, dhdt, dndt] = obj.derivatives(V(i) + k2_V/2, m(i) + k2_m/2, ...
                    h(i) + k2_h/2, n(i) + k2_n/2, time(i) + dt/2, I_stim);
                k3_V = dt * dVdt;
                k3_m = dt * dmdt;
                k3_h = dt * dhdt;
                k3_n = dt * dndt;
                
                [dVdt, dmdt, dhdt, dndt] = obj.derivatives(V(i) + k3_V, m(i) + k3_m, ...
                    h(i) + k3_h, n(i) + k3_n, time(i) + dt, I_stim);
                k4_V = dt * dVdt;
                k4_m = dt * dmdt;
                k4_h = dt * dhdt;
                k4_n = dt * dndt;
                
                V(i+1) = V(i) + (k1_V + 2*k2_V + 2*k3_V + k4_V) / 6;
                m(i+1) = m(i) + (k1_m + 2*k2_m + 2*k3_m + k4_m) / 6;
                h(i+1) = h(i) + (k1_h + 2*k2_h + 2*k3_h + k4_h) / 6;
                n(i+1) = n(i) + (k1_n + 2*k2_n + 2*k3_n + k4_n) / 6;
            end
        end
        
        % 결과 플롯
        function plot_results(obj, time, V, m, h, n)
            figure('Position', [100 100 800 600]);
            
            subplot(2,1,1);
            plot(time, V, 'LineWidth', 2);
            title('Hodgkin-Huxley 모델 시뮬레이션');
            ylabel('막전위 (mV)');
            grid on;
            
            subplot(2,1,2);
            plot(time, m, 'r-', time, h, 'g-', time, n, 'b-', 'LineWidth', 2);
            xlabel('시간 (ms)');
            ylabel('게이팅 변수');
            legend('m (Na+ 활성화)', 'h (Na+ 비활성화)', 'n (K+ 활성화)');
            grid on;
        end
    end
end