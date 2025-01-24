Данный проект содержит код **Scilab** для моделирования движения наночастиц в ограниченном 2D пространстве с учётом столкновений со стенками и распределений скоростей (максвелловского распределения), а также случайным выбором длины свободного пробега. 

Цель программы — сымитировать движение наночастиц в ограниченном объёме (2D-бокс) при заданных физических параметрах:
- Температура (через константу Больцмана k и массу частицы m);
- Начальные координаты и скорости наночастиц;
- Число итераций моделирования;
- Длина свободного пробега, выбираемая случайным образом;
- Столкновений со стенками и время «залипания» на стенке.

  Основой является классическая модель свободного пробега частицы между столкновениями, а скорости задаются через среднеквадратичное (сигму), вытекающее из распределения Максвелла.

  ## Краткое объяснение основных функций

**get_params()**
Функция инициализирует параметры системы и возвращает два списка:
1. **Particle_Par** — параметры частицы.
   - Начальные координаты coord_start.
   - Начальная скорость vel.
   - Диаметр частицы d.
   - Масса частицы m.
   - Температура T.

2. **Modelling_Par** — основные параметры моделирования. 
   - Wall_seg: число делений у стенки (условный параметр для дискретизации)  
   - K: количество итераций основного цикла моделирования  
   - t_lim: лимит по времени (непосредственная остановка моделирования)  
   - H: число частиц (или запусков)  
   - RNG_Av: мат. ожидание распределения (для случайных величин)
   - sigma: стандартное отклонение, связываемое с распределением Максвелла  
   - Lengths: массив длин свободного пробега  
   - TAU: время залипания на стенке  
   - x_lim, y_lim: границы бокса по осям X и Y

Основная функция моделирования **begit_anjumanya(Particle_Par, Modelling_Par)** выполняет:
- Распаковку переданных параметров.
- Основной цикл по количеству шагов K.  
- Вычисление времени между столкновениями и новой координаты частицы.  
- Учёт столкновений со стенками, «залипания» (TAU) и выход за пределы бокса.  
- Сбор результатов (координат, скоростей и пр.).  
- Возвращает структуру или список, содержащий историю движения частиц, результаты столкновений, флаги перелётов и т.д.

Реализованы генераторы случайных величин random_dist_one_num() и random_dist_array() по методу усечения равномерного распределения.
