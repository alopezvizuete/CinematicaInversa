function Ej221_run_IK()
    %Alejandro Lopez Vizuete
    %MIGJRV - Expediente 410
    
    close all;
    %NUMEROS DE TARGET QUE SE MOSTRARAN
    max_target = 5;
    
    % Creamos los elementos del skeleton
    trans0 = [0,0];
    trans1 = [-1,-1]; 
    trans2 = [0,-1];  
    trans3 = [0,-1];  
    
    % Damos valor a las rotaciones iniciales
     rot1 = 0;         
     rot2 = 0;          
     rot3 = 0;  
     
    %Creamos el array skeleton
    skeleton = [trans0; trans1; trans2; trans3];
    
    %Creamos el array rots con las rotaciones iniciales
    rots = [rot1 rot2 rot3];
    
    %inicializamos el numero de elementos que vamos a buscar
    n_elementos = 0;
    
    %Bucle de max_target(especificado en las primeras lineas) elementos a buscar
    while(n_elementos < max_target)
        %Pintamos la grafica
        figure;
        hold on;
        grid on;
        axis equal;
        xlim([-5 5]);
        ylim([-5 5]);

        %Se escogen al azar dos numeros aleatorios
        %que seran el target a tocar
        n1=round(4*rand);
        n2=round(4*rand);
        target = [ 1 0 n1;
        0 1 n2;
        0 0 1];
        %Pintamos el target
        scatter(n1, n2,30,'x');
        %Calculamos la posicion actual
        current_position = compute_position(skeleton, rots);
        %Inicializamos la variable de iteracciones para tener un limite.
        iterac = 0;
        
        %Bucle para aumentar los angulos y llegar al target
        %Se para cuando la distancia entre target y posicion actual es
        %inferior a 0.1, o si llega a 50 iteraciones
        while(distancia(current_position,target)>0.1) && (iterac < 50)
            %Aumentamos el numero de iteraciones
            iterac = iterac +1;
            %Pintamos el target
            scatter(n1, n2,30,'x');
            %Pintamos el skeleton
            plot_skeleton(skeleton, rots);
            drawnow;
            
            %Calculamos la matriz J mediante la tecnica del Jacobiano
            J = compute_jacobian(skeleton,rots);
            %Calculamos su pseudoinversa
            Jinv = pinv(J);
            %Guardamos la posicion actual
            current_position = compute_position(skeleton, rots);
            %Creamos una matriz con la diferencia entre el target y la
            %posicion actual
            diferen = [target(1,3)-current_position(1,3),target(2,3)-current_position(2,3)];
            % Y la multiplicamos por la matriz pseudoinversa jacobiana para
            % obtener los aumentos de angulos
            update = Jinv * transpose(diferen);

            %Una vez tengamos los incrementos de angulos, los sumamos al
            %angulo anterior para poder movernos hacia el target.
            rots(1) = rots(1)+ update(1)*0.2;
            rots(2) = rots(2)+ update(2)*0.2;
            rots(3) = rots(3)+ update(3)*0.2;
            cla;
        end
        %Cuando acabe de todas las iteraciones para el elemento, pasamos al
        %siguiente elemento, volviendo a calcular un target aleatorio
        n_elementos = n_elementos +1;
    end

end

function plot_skeleton(skeleton, rots)
    %Funcion dada en los ejercicios
    rot1 = rots(1);
    rot2 = rots(2);
    rot3 = rots(3);
    % Aplicamos nodo 0

    rot0_matrix = [ 1 0 0;
    0 1 0;
    0 0 1];

    trans0_matrix = [1 0 skeleton(1,1);
    0 1 skeleton(1,2);
    0 0 1];

    pos             = rot0_matrix * trans0_matrix;
    scatter(pos(1,3), pos(2,3),300,'o','filled');

    % Aplicamos nodo 1
    rot1_matrix = [ +cosd(rot1) -sind(rot1) 0;
    +sind(rot1) +cosd(rot1) 0;
    0          0          1 ];

    trans1_matrix = [1 0 skeleton(2,1);
    0 1 skeleton(2,2);
    0 0 1];

    old_pos = pos;
    pos     =   rot0_matrix * trans0_matrix ...
    * rot1_matrix * trans1_matrix ;

    scatter(pos(1,3), pos(2,3),300,'o','filled');
    line([pos(1,3),old_pos(1,3)],[pos(2,3),old_pos(2,3)], 'LineWidth', 2);
    
    % Aplicamos nodo 2
    rot2_matrix = [ +cosd(rot2) -sind(rot2) 0;
    +sind(rot2) +cosd(rot2) 0;
    0          0          1 ];

    trans2_matrix = [1 0 skeleton(2,1);
    0 1 skeleton(2,2);
    0 0 1];

    old_pos = pos;
    pos     =   pos...
    * rot2_matrix * trans2_matrix ;

    scatter(pos(1,3), pos(2,3),300,'o','filled');
    line([pos(1,3),old_pos(1,3)],[pos(2,3),old_pos(2,3)], 'LineWidth', 2);
    
    % Aplicamos nodo 3
    rot3_matrix = [ +cosd(rot3) -sind(rot3) 0;
    +sind(rot3) +cosd(rot3) 0;
    0          0          1 ];

    trans3_matrix = [1 0 skeleton(2,1);
    0 1 skeleton(2,2);
    0 0 1];

    old_pos = pos;
    pos     =   pos...
    * rot3_matrix * trans3_matrix ;

    scatter(pos(1,3), pos(2,3),300,'o','filled');
    line([pos(1,3),old_pos(1,3)],[pos(2,3),old_pos(2,3)], 'LineWidth', 2);
    

end

function pos = compute_position(skeleton, rots)
    %Misma funcion que plot_skeletion, pero devolviendo el valor de 
    % la posicion en el ultimo elemento
    rot1 = rots(1);
    rot2 = rots(2);
    rot3 = rots(3);
    % Aplicamos nodo 0

    rot0_matrix = [ 1 0 0;
    0 1 0;
    0 0 1];

    trans0_matrix = [1 0 skeleton(1,1);
    0 1 skeleton(1,2);
    0 0 1];

    pos             = rot0_matrix * trans0_matrix;

    % Aplicamos nodo 1
    rot1_matrix = [ +cosd(rot1) -sind(rot1) 0;
    +sind(rot1) +cosd(rot1) 0;
    0          0          1 ];

    trans1_matrix = [1 0 skeleton(2,1);
    0 1 skeleton(2,2);
    0 0 1];

    old_pos = pos;
    pos     =   rot0_matrix * trans0_matrix ...
    * rot1_matrix * trans1_matrix ;
   
    % Aplicamos nodo 2
    rot2_matrix = [ +cosd(rot2) -sind(rot2) 0;
    +sind(rot2) +cosd(rot2) 0;
    0          0          1 ];

    trans2_matrix = [1 0 skeleton(2,1);
    0 1 skeleton(2,2);
    0 0 1];

    old_pos = pos;
    pos     =   pos...
    * rot2_matrix * trans2_matrix ;
    
    % Aplicamos nodo 3
    rot3_matrix = [ +cosd(rot3) -sind(rot3) 0;
    +sind(rot3) +cosd(rot3) 0;
    0          0          1 ];

    trans3_matrix = [1 0 skeleton(2,1);
    0 1 skeleton(2,2);
    0 0 1];

    old_pos = pos;
    pos     =   pos...
    * rot3_matrix * trans3_matrix ;
end

function dist = distancia(current_position,target)
    %Funcion para calcular la distancia entre dos puntos dados en matrices
    dist = sqrt((current_position(1,3)-target(1,3)).^2+(current_position(2,3)-target(2,3)).^2);
end

function j = compute_jacobian(skeleton,rots)
    %Funcion para obtener la matriz del Jacobino
    %Obtenemos la posicion inicial del ultimo elemento al entrar a la funcion
    pos_inicial = compute_position(skeleton, rots);
    
    %Aumentamos el valor de la primera rotacion en 2 unidades (atenuandola
    %al calcular el update, ya que deben ser cambios pequeños de rotacion)
    rots(1) = rots(1) + 2;
    %Obtenemos la posicion del ultimo elemento al aumentar esta rotacion
    pos_sumAng = compute_position(skeleton, rots);
    %Creamos una matriz que contenga la diferencia entre la posicion con el
    %desplazamiento y la posicion inicial (incremento de distancia)
    pos_res =  [ 1 0 pos_sumAng(1,3)-pos_inicial(1,3);
                0 1 pos_sumAng(2,3)-pos_inicial(2,3);
                0 0 1];
    %El incremento de distancia lo introducimos en la matriz J, en su
    %espacio correspondiente dividido entre el angulo que hemos sumado
    %anteriormente (incremento de angulo)
    j= [ pos_res(1,3)/2 0 0;
       pos_res(2,3)/2 0 0];
   
    %Y dejamos el angulo como estaba
    rots(1) = rots(1) - 2;
    
    %Este paso lo repetimos para el angulo 2 y 3
    rots(2) = rots(2) + 2;
    pos_sumAng = compute_position(skeleton, rots);
    pos_res =  [ 1 0 pos_sumAng(1,3)-pos_inicial(1,3);
                0 1 pos_sumAng(2,3)-pos_inicial(2,3);
                0 0 1];
    j(1,2) = pos_res(1,3)/2;
    j(2,2) = pos_res(2,3)/2;
   
    rots(2) = rots(2) - 2;
    
    rots(3) = rots(3) + 2;
    pos_sumAng = compute_position(skeleton, rots);
    pos_res =  [ 1 0 pos_sumAng(1,3)-pos_inicial(1,3);
                0 1 pos_sumAng(2,3)-pos_inicial(2,3);
                0 0 1];
    j(1,3) = pos_res(1,3)/2;
    j(2,3) = pos_res(2,3)/2;
    rots(3) = rots(3) - 2;
end