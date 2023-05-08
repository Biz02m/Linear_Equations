clear all
close all
%188594
c = 9;
d = 4;
e = 5;
N = 9*100 + c*10 + d; 
a1 = 5 + e;
f = 8;

%Zadanie A
A1 = DiagTrian(scalarMult(ones(N),a1));
A2 = DiagTrian(scalarMult(ones(N), -1), 1);
A3 = DiagTrian(scalarMult(ones(N), -1), 2);


A = A1 + A2 + A3 + MyTransposition(A2) + MyTransposition(A3); %utworzenie macierzy A

b = ones(N,1);

for i = 1:N
    b(i) = sin(i*(f+1)); %utworzenie wektora b
end

%-----------------------------
%Zadanie B
disp("Zadanie B")

%Metoda Jacobiego
[Jacobi_iterations, Jacobi_czas, Jacobi_norm, Jacobi_norms] = Jacobi(A,b,N);

fprintf("Liczba Iteracji Metodą Jacobiego\n");
disp(Jacobi_iterations);
fprintf("Czas wykonywania Metodą Jacobiego\n");
disp(Jacobi_czas);
fprintf("Norma błędu residuum Metodą Jacobiego\n");
disp(Jacobi_norm);
figure
semilogy(Jacobi_norms), title("Wykres norm błędu residuum w iteracjach Jacobiego Zadanie B"), xlabel("Liczba iteracji"), ylabel("Norma błędu residuum"), legend("Metoda Jacobiego");

fprintf("-----------------\n");
%Metoda Gaussa Siedla
[GS_iterations, GS_czas, GS_norm, GS_norms] = Gauss_Seidl(A,b,N);

fprintf("Liczba Iteracji Metodą Gaussa Seidla\n");
disp(GS_iterations);
fprintf("Czas wykonywania Metodą Gaussa Seidla\n");
disp(GS_czas);
fprintf("Norma błędu residuum Metodą Gaussa Seidla\n");
disp(GS_norm);
figure
semilogy(GS_norms), title("Wykres norm błędu residuum w iteracjach Gaussa-Seidla Zadanie B"), xlabel("Liczba iteracji"), ylabel("Norma błędu residuum"), legend("Gauss-Seidl");



%---------------------------------
%zadanie C i D
disp("========================")
disp("zadanie C i D")

clear all

c = 9;
d = 4;
N = 9*100 + c*10 + d;
f = 8;

A1 = DiagTrian(scalarMult(ones(N),3));
A2 = DiagTrian(scalarMult(ones(N), -1), 1);
A3 = DiagTrian(scalarMult(ones(N), -1), 2);

A = A1 + A2 + A3 + MyTransposition(A2) + MyTransposition(A3); %utworzenie macierzy A

b = ones(N,1);

for i = 1:N
    b(i) = sin(i*(f+1)); %utworzenie wektora b
end


%Metoda Bezpośrednia Faktoryzacji LU
[LU_norm, LU_time] = ImmediateLU(A,b);

fprintf("Liczba Iteracji Metodą Faktoryzacji LU\n");
disp("1 ponieważ ta metoda jest bezpośrednia :)");
fprintf("Czas wykonywania Metodą Faktoryzacji LU\n");
disp(LU_time);
fprintf("Norma błędu residuum Metodą Faktoryzacji LU\n");
disp(LU_norm);

fprintf("-----------------\n");
%Metoda Jacobiego
[Jacobi_iterations, Jacobi_czas, Jacobi_norm, Jacobi_norms] = Jacobi(A,b,N);

fprintf("Liczba Iteracji Metodą Jacobiego\n");
disp(Jacobi_iterations);
fprintf("Czas wykonywania Metodą Jacobiego\n");
disp(Jacobi_czas);
fprintf("Norma błędu residuum Metodą Jacobiego\n");
disp(Jacobi_norm);
figure 
semilogy(Jacobi_norms), title("Wykres norm błędu residuum w iteracjach Jacobiego Zadanie C"), xlabel("Liczba iteracji"), ylabel("Norma błędu residuum"), legend("Metoda Jacobiego");


fprintf("-----------------\n");
%Metoda Gaussa Siedla
[GS_iterations, GS_czas, GS_norm, GS_norms] = Gauss_Seidl(A,b,N);

fprintf("Liczba Iteracji Metodą Gaussa Seidla\n");
disp(GS_iterations);
fprintf("Czas wykonywania Metodą Gaussa Seidla\n");
disp(GS_czas);
fprintf("Norma błędu residuum Metodą Gaussa Seidla\n");
disp(GS_norm);
figure 
semilogy(GS_norms), title("Wykres norm błędu residuum w iteracjach Gauss-Seidl Zadanie C"), xlabel("Liczba iteracji"), ylabel("Norma błędu residuum"), legend("Gauss-Seidl");

%Zadanie E

disp("========================")
disp("zadanie E")
clear all

c = 9;
d = 4;
N = [500, 1000, 1500 , 2000, 2500, 3000, 3500];
f = 8;
c = 9;
d = 4;
e = 5;
a1 = 5 + e;
f = 8;

Jacobi_czasy = zeros(1, numel(N));
GS_czasy = zeros(1, numel(N));
LU_czasy = zeros(1, numel(N));

for i = 1:numel(N)
    A1 = DiagTrian(scalarMult(ones(N(i)),a1));
    A2 = DiagTrian(scalarMult(ones(N(i)), -1), 1);
    A3 = DiagTrian(scalarMult(ones(N(i)), -1), 2);
    A = A1 + A2 + A3 + MyTransposition(A2) + MyTransposition(A3); %utworzenie macierzy A

    b = ones(N(i),1);

    for j = 1:N(i)
        b(j) = sin(j*(f+1)); %utworzenie wektora b
    end

    [LU_norm, LU_time] = ImmediateLU(A,b);
    LU_czasy(i) = LU_time;
    [Jacobi_iterations, Jacobi_czas, Jacobi_norm, Jacobi_norms] = Jacobi(A,b,N(i));
    Jacobi_czasy(i) = Jacobi_czas;
    [GS_iterations, GS_czas, GS_norm, GS_norms] = Gauss_Seidl(A,b,N(i));
    GS_czasy(i) = GS_czas;
    disp(i)
end

figure
plot(N,LU_czasy);
hold on
plot(N,Jacobi_czasy, "b- ");
plot(N,GS_czasy, "r- "), xlabel("Liczba niewiadomych"),ylabel("Czas Obliczeń [s]"), title("Czas rozwiązywania układów równań metodami iteracyjnymi"), legend("Faktoryzacja LU", "Jacobi", "Gauss-Seidl");
