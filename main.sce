clear; xdel(winsid())
///////////////////////////////////////////////////////////////////////////////
//--------------------------ИСХОДНЫЕ ДАННЫЕ----------------------------------//
///////////////////////////////////////////////////////////////////////////////
function [Particle_Par, Modelling_Par] = get_params() //входные параметры основного цикла
    //константы
    k = 1.38 * 10^-23//константа боцмана
    
    //параметры частицы
    coord_start = [random_dist_one_num(0,10,10) random_dist_one_num(0,10,10)] //начальные координаты
    vel = [0.000001 0.000001] //скорость
    d = 1 //диаметр частицы
    m = 0.0000000001 //масса частицы, г
    T = 10 //температура, град К
    
    //настройки
    Wall_seg = 50 //число делений у стенки
    K = 1000 // количество циклов
    t_lim = 5 * 10^8 //предел 
    H = 2 // число частиц
    RNG_Av = 0 // мат ожидание распределения
    sigma = sqrt(k*T/m) //распределение максвелла
    Lengths = random_dist_array(K, 0, 10, 10) //длина свободного пробега
    TAU = 10^6 //время зависания на стенке
    //границы (стенки)
    x_lim = [0 100]
    y_lim = [0 100]
    
    //формирование листов
    Particle_Par = list(coord_start, vel, d, m, T)
    Modelling_Par = list(Wall_seg, K, t_lim, H, RNG_Av, sigma, Lengths, TAU, x_lim, y_lim)
endfunction

///////////////////////////////////////////////////////////////////////////////
//------------------------------- ФУНКЦИИ -----------------------------------//
///////////////////////////////////////////////////////////////////////////////
//основной цикл

function Model_output = begit_anjumanya(Particle_Par, Modelling_Par)
    //распаковка листов
    //параметры частицы
    coord = Particle_Par(1) //начальные координаты
    vel = Particle_Par(2) //скорость
    d = Particle_Par(3) //диаметр частицы
    m = Particle_Par(4) //масса частицы, г
    T = Particle_Par(5) //температура, град К
    
    //параметры моделирования
    Wall_seg = Modelling_Par(1) //число делений у стенки
    K = Modelling_Par(2) // количество циклов
    t_lim = Modelling_Par(3) //предел 
    H = Modelling_Par(4) // число частиц
    RNG_Av = Modelling_Par(5) // мат ожидание распределения
    sigma = Modelling_Par(6) //распределение максвелла
    Lengths = Modelling_Par(7) //длина свободного пробега
    TAU = Modelling_Par(8) //время зависания на стенке
    x_lim = Modelling_Par(9) //границы (стенки)
    y_lim = Modelling_Par(10)
    
    //пустые векторы для результатов
    v_up = [] 
    v_down = []
    v_left = []
    v_right = []
    x_vrezalsya_up = []
    x_vrezalsya_down = []
    y_vrezalsya_left = []
    y_vrezalsya_right = []
    x_perelet_flag = %F
    y_perelet_flag = %F
    t_max = 0
    coord_new = [0 0] //для новых координат
    coord_stop_x = [0 0]
    coord_stop_y = [0 0]
    coord_arch = [0 0] //архив координат

    for n = 1:1:K
        vel_sum = sqrt(sum(vel.^2)) //суммарная скорость
        L = Lengths(n) //длина свободного пробега
        t = L/vel_sum //время, нс
        coord_prev = coord //предыдущая координата
        coord = coord($,:) //берем последнюю строку
        coord = coord + vel.*t  //изменение координат
        
//-----------------------------------------------------------------------------
//УЧЕТ ГРАНИЦ -----------------------------------------------------------------
//-----------------------------------------------------------------------------
        //для Х 
        if coord(1) < x_lim(1) then
            x_perelet_flag = %T //флаг перелета
            perelet = abs(x_lim(1) - coord(1))
            x_vrezalsya = x_lim(1) //точка x в которой врезался
            y_vrezalsya_left = get_y(x_vrezalsya, coord(2), coord(1), coord_prev(2), coord_prev(1)) //точка y в которой врезался
            reflected_x_left =  coord(1) + 2*perelet //оценочно куда модекула отлетит после отражения, бывш coord_new(1)
            
            //зависание частицы
            t = t + TAU //учет времени зависания
            coord_stop_x(1) = x_lim(1) //координаты в которых зависла частица
            coord_stop_x(2) = y_vrezalsya_left
            
            //координата после отражения
            vel = rand_new_vel(Modelling_Par) //новая скорость
            while vel(1) < 0
                vel = rand_new_vel(Modelling_Par)
            end
            coord_new(1) = coord_stop_x(1) + vel(1).*t
        end
    
        if coord(1) > x_lim(2) then
            x_perelet_flag = %T
            perelet = abs(x_lim(2) - coord(1))
            x_vrezalsya = x_lim(2)
            y_vrezalsya_right = get_y(x_vrezalsya, coord(2), coord(1), coord_prev(2), coord_prev(1))
            reflected_x_right = coord(1) - 2*perelet
            
            //зависание частицы
            t = t + TAU //учет времени зависания
            coord_stop_x(1) = x_lim(2) //координаты в которых зависла частица
            coord_stop_x(2) = y_vrezalsya_right
            
            //координата после отражения
            vel = rand_new_vel(Modelling_Par) //новая скорость
            while vel(1) > 0
                vel = rand_new_vel(Modelling_Par)
            end
            coord_new(1) = coord_stop_x(1) + vel(1).*t
        end 
    
        //для Y 
        if coord(2) < y_lim(1) then
            y_perelet_flag = %T
            perelet = abs(y_lim(1) - coord(2))
            y_vrezalsya = y_lim(1)
            x_vrezalsya_down = get_x(y_vrezalsya, coord(2), coord(1), coord_prev(2), coord_prev(1))
            reflected_y_down =  coord(2) + 2*perelet
            
            //зависание частицы
            t = t + TAU //учет времени зависания
            coord_stop_y(1) = y_lim(1) //координаты в которых зависла частица
            coord_stop_y(2) = x_vrezalsya_down
            
            //координата после отражения
            vel = rand_new_vel(Modelling_Par) //новая скорость
            while vel(2) < 0
                vel = rand_new_vel(Modelling_Par)
            end
            coord_new(2) = coord_stop_y(2) + vel(2).*t
        end 
        
        if coord(2) > y_lim(2) then
            y_perelet_flag = %T
            perelet = abs(y_lim(2) - coord(2))
            y_vrezalsya = y_lim(2)
            x_vrezalsya_up = get_x(y_vrezalsya, coord(2), coord(1), coord_prev(2), coord_prev(1))
            reflected_y_up =  coord(2) - 2*perelet
            
            //зависание частицы
            t = t + TAU //учет времени зависания
            coord_stop_y(1) = y_lim(2) //координаты в которых зависла частица
            coord_stop_y(2) = x_vrezalsya_up

            //координата после отражения
            vel = rand_new_vel(Modelling_Par) //новая скорость
            while vel(2) > 0
                vel = rand_new_vel(Modelling_Par)
            end
            coord_new(2) = coord_stop_y(2) + vel(2).*t
        end 
        //учитываем
        if (x_perelet_flag && y_perelet_flag) == %T then
            coord(1) = coord_new(1)
            coord(2) = coord_new(2)
            coord = cat(1, coord_stop_x, coord_stop_y, coord)
        end
        if x_perelet_flag == %T then
            coord(1) = coord_new(1)
            coord = cat(1, coord_stop_x, coord)
        end
        if y_perelet_flag == %T then
            coord(2) = coord_new(2)
            coord = cat(1, coord_stop_y, coord)
        end
        disp("Цикл")
        disp(n)
        disp(coord)
        
//-----------------------------------------------------------------------------
//Учет границ -----------------------------------------------------------------
//-----------------------------------------------------------------------------
        
        vel = rand_new_vel(Modelling_Par) //новая скорость
        
        //формируем архивы
        t_arch(1, n) = t //архив времени
        coord_arch = cat(1, coord_arch, coord)
//        coord_arch(1, n) = coord(1) // архив координат 
//        coord_arch(2, n) = coord(2)
        vel_arch(1, n) = vel(1) // архив скоростей
        vel_arch(2, n) = vel(2)
        
        //складируем то что врезалось для построения гистограмм
        v_up = cat(1, v_up, x_vrezalsya_up) 
        v_down = cat(1, v_down, x_vrezalsya_down) 
        v_left = cat(1, v_left, y_vrezalsya_left) 
        v_right = cat(1, v_right, y_vrezalsya_right) 
        
        //зануляем чтобы не перескладировалос
        x_vrezalsya_up = []
        x_vrezalsya_down = []
        y_vrezalsya_left = []
        y_vrezalsya_right = []
        coord_new = []
        x_perelet_flag = %F
        y_perelet_flag = %F
        
        //остановить при достижении заданного времени
        t_max = t_max + t
        if t_max > t_lim
            n = K + 1
        end
    end
    
    //формирование выходного листа
    Model_output = list(coord_arch, v_up, v_down, v_left, v_right, t_arch, vel_arch)
endfunction

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//ФОРМИРОВАНИЕ РАСПРЕДЕЛЕНИЯ НУЖНОЙ ФОРМЫ
// функция распределения
function y = distribution_func(lambda, L)
    // формула распределения
    y = ((4.*lambda.^2)./(sqrt(%pi).*L^3)).*exp(-(lambda./L).^2)
endfunction

// формирование одного числа
function num = random_dist_one_num(x_min, x_max, L)
    temp = 1
    while temp < 2 
        x = grand(1, 1, "unf", x_min, x_max)
        y = grand(1, 1, "unf", x_min, x_max)
        if y < distribution_func(x, L) then
            num = x
            temp = 4 //это я так цикл останавливаю
        end
    end
endfunction

// формирование массива чисел
function array = random_dist_array(num_of_elements, x_min, x_max, L)
    //num_of_elements - число элементов в массиве
    //x_min, x_max - границы усечения, напр 0 и 10
    //L - параметр распределения, напр 10
    i = 1
    while i <= num_of_elements
        array(:,i) = random_dist_one_num(x_min, x_max, L)
        i = i + 1
    end
endfunction

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//расчет новых координат

function x = get_x(y,y1,x1,y2,x2)
    //y1,y2 - имеющиеся координаты
    //x1,x2 - тоже имеющиеся
    //y - координата, для которой рассчитывается x
    x = ((y - y1)/(y2 - y1))*(x2 - x1) + x1
endfunction

function y = get_y(x,y1,x1,y2,x2)
    //y1,y2 - имеющиеся координаты
    //x1,x2 - тоже имеющиеся
    //x - координата, для которой рассчитывается y
    y = ((x - x1)/(x2 - x1))*(y2 - y1) + y1
endfunction

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//рандом новой скорости 

function velocity = rand_new_vel(Modelling_Par)
    RNG_Av = Modelling_Par(5)
    sigma = Modelling_Par(6)
    velocity = grand(1, 2, "nor", RNG_Av, sigma) //новая скорость
endfunction

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//histc для версии 6.0.2
//зачем вообще аргументы местами менять?????

function histcout = histc_sciver(arg1, arg2)
    scilab_ver = getversion()
    if scilab_ver == "scilab-6.1.1" then
        histcout = histc(arg1, arg2)
    else
        histcout = histc(arg2, arg1)
    end
endfunction

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//--------------------------------- НАЧАЛО ----------------------------------//
///////////////////////////////////////////////////////////////////////////////

//ИНИЦИАЛИЗАЦИЯ
[Particle_Par, Modelling_Par] = get_params() //получаем исходные данные
H = Modelling_Par(4) // число частиц
Wall_seg = Modelling_Par(1) //число делений у стенки
v_up_all = []
v_down_all = []
v_left_all = [] 
v_right_all = []

//МОДЕЛИРОВАНИЕ
f = figure(1) //создаем окно
for nn = 1:1:H 
    Model_output = begit_anjumanya(Particle_Par, Modelling_Par) //запуск моделирования
    //распаковка листа
    coord_arch = Model_output(1)'
    v_up = Model_output(2)
    v_down = Model_output(3)
    v_left = Model_output(4)
    v_right = Model_output(5)
    
    //рисуем графики
    comet(coord_arch(1,:), coord_arch(2,:))
    plot(coord_arch(1,:), coord_arch(2,:))
    last_line = gce() //настраиваем цвета графиков
    last_line.children.foreground = nn
    
    //считаем общее число ударов
    v_up_all = cat(1, v_up_all, v_up) 
    v_down_all = cat(1, v_down_all, v_down) 
    v_left_all = cat(1, v_left_all, v_left) 
    v_right_all = cat(1, v_right_all, v_right) 
end
gca().grid=[1 1 1]

//строим гистограммы
//верхняя
f_up = figure(2)
title("Верхняя стенка")
up_hist = histc_sciver(v_up_all, Wall_seg) 
bar(up_hist)
//нижняя
f_down = figure(3)
title("Нижняя стенка")
down_hist = histc_sciver(v_down_all, Wall_seg) 
bar(down_hist)
//левая
f_left = figure(4)
title("Левая стенка")
left_hist = histc_sciver(v_left_all, Wall_seg) 
bar(left_hist)
//правая
f_right = figure(5)
title("Правая стенка")
right_hist = histc_sciver(v_right_all, Wall_seg) 
bar(right_hist)
