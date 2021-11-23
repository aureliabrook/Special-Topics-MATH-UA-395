% Define the Pauli matrices
X = [0 1; 1 0]; % Quantum equivalent of a NOT gate 
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];

% Other quantum logic gates
I = [1 0; 0 1]; % Identity matrix
H = (1/sqrt(2) * (X + Z)); % Hadamard Gate
CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]; % Controlled NOT gate

C = [0 1; 1 0]; %Correction gate

y = makeGIF({{X, H, I}, {X, H, C}}, {[1, 0], [1, 0]}, 20, ["No noise", "Noisy"], "WithCorrection", [false, true]);

function P = P(p) 
    P = [1 0; 0 cos(p)+1i*sin(p)]; 
end % Phase shift with regards to some angle phi

function [euc, ang] = QubitToEuclidean(q)
    a = real(q(1));
    b = imag(q(1));
    if a == 0
        gamma = 0;
    else
        gamma = atan(b / a); 
    end
    theta = 2*acos(a / cos(gamma));
    if theta == 0
        phi = 0;
    else
        phi = acos(real(q(2)) / sin(theta/2)) - gamma;
    end
    
    % Convert Spherical to Euclidean    
    euc = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
    ang = [theta, phi];
end

function q1 = noise(q0)
     n = randi(6);
     if n < 4
         q1 = [1, 0];
         disp("ket0");
     elseif n < 6
        q1 = [1/sqrt(2), 1/sqrt(2)];
        disp("ket+");
      else
          q1 = q0;
          disp("nothing");
     end
end

function b = measure(euc)
    b = [0, 0];
    if euc(3) > 0
        b(3) = 1;
    else
        b(3) = -1;
    end
end

function  q1 = makeGIF(gates, q0, n, ti, filename, err)
    q1 = {};
    vecs = {{}};
    d = length(q0);
    for i = 1:d
        for k = 1:length(gates{i})
            [euc0, ~] = QubitToEuclidean(q0{i});
            q1{i} = gates{i}{k} * q0{i}';
            [euc1, ~] = QubitToEuclidean(q1{i});
            if k == 1
                vecs{i} = {euc0};
            else 
                vecs{i}{end+1} = euc0;
            end
            vec_diff = (euc1 - euc0) / n;
            if round(euc0*euc1') == -1
                for j = 1:n-1
                    vecs{i}{end+1} = (euc0+j*vec_diff);
                end
            else
                for j = 1:n-1
                    vecs{i}{end+1} = (euc0+j*vec_diff)/norm(euc0+j*vec_diff);
                end
            end
            vecs{i}{end+1}=euc1;
            q0{i} = q1{i}';
            if err(i)
                q1{i} = noise(q0{i})';
                euc0 = euc1;
                [euc1, ~] = QubitToEuclidean(q1{i});
                vecs{i}{end+1} = euc0;
                vec_diff = (euc1 - euc0) / n;
                if round(euc0*euc1') == -1
                    for j = 1:n-1
                        vecs{i}{end+1} = (euc0+j*vec_diff);
                    end
                else
                    for j = 1:n-1
                        vecs{i}{end+1} = (euc0+j*vec_diff)/norm(euc0+j*vec_diff);
                    end
                end
                vecs{i}{end+1}=euc1;
                q0{i} = q1{i}';
            elseif ~err(i) && ~all(~err)
                for j = 1:n+1
                    vecs{i}{end+1} = euc1;
                end
            end
        end
    end
    for i = 1:d
        vecs{i}{end+1} = vecs{i}{end};
        euc0 = vecs{i}{end};
        euc1 = measure(euc0);
        vec_diff = (euc1 - euc0) / n;
        if round(euc0*euc1') == -1
            for j = 1:n-1
                vecs{i}{end+1} = (euc0+j*vec_diff);
            end
        else
            for j = 1:n-1
                vecs{i}{end+1} = (euc0+j*vec_diff)/norm(euc0+j*vec_diff);
            end
        end
        vecs{i}{end+1}=euc1;
    end
    [x, y, z] = sphere(25);
    h = figure('visible', 'off');
    for i = 1:length(vecs{1})
        % Create the Bloch sphere, and customize how it looks
        for j = 1:d
            subplot(1, d, j);
            bloch = surf(x, y, z);
            title(ti(j));
            axis equal manual
            colormap autumn
            shading interp
            bloch.FaceAlpha = 0.1;
            if ~err(j) && rem(i, 2*n+2) > n+1 && ~all(~err)
                color = "Green";
            elseif 1 < i && i < length(vecs{1}) && rem(i, n+1) ~= 0
                color = "Red";
            else
                color = "Blue";
            end
            if (i > 2 * length(gates{d}) * (n + 1) && ~all(~err)) || (i > length(gates{d}) * (n + 1) && all(~err))
                color = "Black";
            end
                
            % Create axes and labels
            line([-1, 1], [0, 0], [0, 0], 'LineWidth', 1, 'Color', [0 0 0]);
            line([0, 0], [-1, 1], [0, 0], 'LineWidth', 1, 'Color', [0 0 0]);
            line([0, 0], [0, 0], [-1, 1], 'LineWidth', 1, 'Color', [0 0 0]);
            text(0, 0, 1.25, '$\left| 0 \right>$', 'Interpreter', 'latex');
            text(0, 0, -1.25, '$\left| 1 \right>$', 'Interpreter', 'latex');
            text(1.25, 0, 0, '$\left| + \right>$', 'Interpreter', 'latex');
            text(-1.25, 0, 0, '$\left| - \right>$', 'Interpreter', 'latex');
            hold on
            quiver3(0, 0, 0, vecs{j}{i}(1), vecs{j}{i}(2), vecs{j}{i}(3), 'LineWidth', 2, 'Color', color);
            hold off
        end
        drawnow;
        
        % Capture the plot as an image 
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 

        % Write to the GIF File 
        if i == 1 
            imwrite(imind,cm,strcat(filename, ".gif"),'gif', 'DelayTime',0.025, 'Loopcount',inf); 
        else 
            imwrite(imind,cm,strcat(filename, ".gif"),'gif','DelayTime',0.025, 'WriteMode','append'); 
        end

        if i == 1 || i == length(vecs{1}) || rem(i, n+1) == 0
            for k = 1:9
                imwrite(imind,cm,strcat(filename, ".gif"),'gif','DelayTime',0.025, 'WriteMode','append'); 
            end
        end
    end
end