# Symulacja paczki falowej napotykającej próg potencjału

Autor: Antonina Pater
<br/><br/>
1. Cel projektu <br/>
Napisanie programu rozwiązującego równanie Schroedingera dla elektronu trafiającego na próg potencjału i porównanie czasu jego wykonania przy jednym i wielu wątkach.

2. Algorytm działania <br/>
W programie rozwiązywane jest zależne od czasu równanie Schroedingera (TDSE) dla gaussowskiej paczki falowej opisującej elektron [2].
Do przeprowadzenia obliczeń zastosowano metodę Cranka-Nicolsona, oryginalnie zastosowaną do rozwiązywania równań różniczkowych cząstkowych dla przewodnictwa ciepła [1].
Dzięki niej otrzymujemy układ równań liniowych: <br/>
A &Psi;(t+dx) = B&Psi;(t), <br/>
gdzie: <br/>
A, B - odpowiednie macierze związane z samym równaniem (w naszym przypadku są one macierzami rzadkimi, zawierają niezerowe elementy tylko na diagonali, nad i pod nią), <br/>
&Psi;(t) - wektor zawierający wartości funkcji falowej w czasie t. <br/>
Takie równanie jest rozwiązywane z pomocą faktoryzacji LU (metoda ta jest szybsza dla macierzy rzadkich). <br/> <br/>
Schemat działania: <br/>
a) inicjalizacja funkcji falowej, <br/>
b) uzupełnienie macierzy A i B oraz dekompozycja macierzy A metodą LU, <br/>
c) rozwiązanie równania, <br/>
d) normalizacja funkcji falowej, <br/>
e) powrót do ppkt. c lub zatrzymanie symulacji. <br/>

3. Budowa programu <br/>
Program składa się z następujących plików:
	* main.py - zawierającego wywołanie menu użytkownika oraz przeprowadzającego obliczenia w pojedynczym wątku oraz w systemie wielowątkowym,
	* sim_params.py - zawierającego zmienne z parametrami symulacji oraz funkcje do obsługi menu,
	* waveeq_cn.py - zawierającego funkcje do przeprowadzenia symulacji, liczenia obserwabli (energia kinetyczna i potencjalna) i normy funkcji falowej, <br/> 
oraz oddzielnego folderu z testami. <br/> 

    By użyć program należy go pobrać z tej strony i pobrać następujące biblioteki (jeśli nie są już zainstlowane): <br/>
numpy, scipy, matplotlib, time, multiprocessing.

4. Przykład działania programu.  <br/>
Program symuluje cząstkę o masie elektronu znajdującą się w układzie o długości 40 nm, w którym próg potencjału występuje w połowie układu oraz ma określoną długość i wysokość.
Użytkownik może ustawić ww. dwa parametry progu, parametry paczki gaussowskiej - długość fali związanej z elektronem, jego położenie początkowe oraz szerokość paczki (sigma, odchylenie standardowe).  <br/>
UWAGA 1. W menu w nawiasach kwadratowych podane są wartości domyślne. W celu ich zaakceptowania wystarczy nacisnąć klawisz ENTER. <br/>
UWAGA 2. Standardowo funkcja falowa jest kreślona dla czasu 0 i 90 fs. Użytkownik może dodać własne wartości.<br/>
UWAGA 3. W celu prawidłowego funkcjonowania symulacji zaleca się stosowanie wysokości progu zbliżonej do eneergii początkowej cząstki lub mniejszej. <br/> <br/>
Wartości domyślne:
    * położenie początkowe cząstki: x<sub>0</sub> =10 nm
    * długość fali: &Lambda; = 5 nm
    * szerokość paczki &sigma; = 2.5 nm
    * wysokość progu V<sub>0</sub> = 0.15 eV
    * długość progu l = 20 nm <br/>
    
- Użycie z wszystkimi dostępnymi wartościami domyślnymi:     <br/>
Po wyświetleniu komunikatu początkowego:     <br/>
Wykonać symulację z parametrami domyślnymi (tak/nie)? [nie]     <br/>
Wpisać "tak" i wcisnąć klawisz ENTER     <br/>
Symulacja wykona się z wartościami domyślnymi.     <br/>

- Swobodna propagacja cząstki: <br/>
Po wyświetleniu komunikatu początkowego: <br/>
Wykonać symulację z parametrami domyślnymi (tak/nie)? [nie] <br/>
Wcisnąć klawisz ENTER. <br/>
Wybrać wartości domyślne położenia początkowego, długości fali i szerokości paczki, wciskając klawisz ENTER. <br/>
Wpisać 0 jako wysokość progu i potwierdzić klawiszem ENTER. <br/>
Długość progu potwierdzić klawiszem ENTER. <br/>
Wybrać wartości domyślne czasu wykreślania wciskająć klawisz ENTER. <br/>
Symulacja wykona się z zadanymi parametrami. <br/> <br/>
##### Źródła <br/>
[1] J. Crank and P. Nicolson, "A practical method for numerical evaluation of solutions of partial differential equations of the heat-conduction type," Mathematical Proceedings of the Cambridge Philosophical Society, vol. 43, no. 1, pp. 50–67, 1947. <br/>
[2] D. M. Sullivan, "Quantum Mechanics for Electrical Engineers," 2012.
