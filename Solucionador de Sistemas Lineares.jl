# -------- IMPORTAÇÕES DE BIBLIOTECAS --------
# A FUNÇÃO lu() FOI FEITA CONSIDERANDO RESTRIÇÃO DE QUE NÃO DEVE EXISTIR ZEROS NA DIAGONAL PRINCIPAL DA MATRIZ U DEVIDO O FATO DO OBJETIVO DO PROJETO DE BUSCAR RESOLVER SISTEMAS LINEARES.
# SUJEITO A MUDANÇAS.
using LinearAlgebra # para importar os comandos, funções de algebra linear.
using Random

# ------------- ALGORITMO PARA RESOLVER SISTEMAS TRIANGULARES INFERIORES -------------------
"""
sti(A, b). ONDE "A" É UMA MATRIZ TRIANGULAR INFERIOR DOS COEFICIENTES DE ORDEM n x n E "b" É O VETOR DOS TERMOS INDEPENDENTES DE ORDEM n x 1.

ESSA FUNÇÃO VAI RESOVER UM SISTEMA LINEAR Ax = b USANDO UM SISTEMA TRIANGULAR INFERIOR.

Ex: A = [1 0 0; 2 3 0; 1 4 5] e b = [1, 3, 6]

sti(A, b) = [1, 0.3333333..., 0.73333333...]

PARA TIRAR A PROVA REAL:

julia) x = sti(A, b)

julia) A*x

O RESULTADO SERÁ O VETOR b.
"""
function sti(matriz, vetor) # o algoritmo de fato que resolve o sistema triangular inferior pelo algoritmo descoberto no livro: "Cálculo Numérico: aprendizagem com Apoio de Software"

    # Toda vez que iniciarmos essa função definimos o vetor solução e o tamanho da matriz.

    tamanho = size(matriz, 1) # fornece a o n° de linhas da matriz, como é quadrada o numero de colunas é igual.

    vetor_sol = zeros(Float64, tamanho) # cria um vetor cheio de zeros com tamanho igual ao número de linhas da matriz que se esta resolvendo. Esse é o vetor que vai armazenar a solução.

    for i in 1:tamanho # fazemos um laço for que depende do tamanho da matriz.

        somatorio = 0.0 # para cada iteração, a soma linear muda.

        for j in 1:i-1 # é o somatório de j indo de 1 até i-1 somando a_ij*x_j
            somatorio += matriz[i, j] * vetor_sol[j]
        end

        # aqui calculamos o valor de x após ter realizado o somatório e armazenamos no vetor solução.
        vetor_sol[i] = (vetor[i] - somatorio) / matriz[i, i]

    end

    return vetor_sol

end

# --------------- FIM DA FUNÇÃO ---------------


# ------------------- ALGORITMO PARA RESOLVER SISTEMAS TRIANGULARES SUPERIORES -------------------------

"""
sts(A, b). ONDE "A" É UMA MATRIZ TRIANGULAR SUPERIOR DOS COEFICIENTES DE ORDEM n x n E "b" É O VETOR DOS TERMOS INDEPENDENTES.

ESSA FUNÇÃO VAI RESOLVER O SISTEMA LINEAR Ax = b VIA SISTEMA TRIANGULAR SUPERIOR.

Ex: A = [1 3 4; 0 4 3; 0 0 -1] e b = [3, -2, 1]

sts(A, b) = [6.25, 0.25, -1.0]

PARA TIRAR A PROVA REAL:

julia) x = sts(A, b)

julia) A*x

O RESULTADO SERÁ O VETOR b.
"""
function sts(matriz, vetor) # o algoritmo de fato que resolve o sistema triangular superior pelo algoritmo descoberto no livro: "Cálculo Numérico: aprendizagem com Apoio de Software"

    # toda vez que iniciamos essa função definimos a variável "tamanho" que recebe o tamanho da matriz e a variável "vetor_sol" que recebe o nosso vetor x solução.

    tamanho = size(matriz, 1) # o size retorna quantas linhas tem essa matriz, se colocar o 2 ele retorna a quantidade de colunas.

    vetor_sol = zeros(Float64, tamanho) # ele é composto com entradas de zeros que são do tipo float64.

    for i in tamanho:-1:1 # laço for que vai diminuindo desde a n-ésima linha até a primeira linha.

        somatorio = 0.0 # a soma linear que

            for j in i+1:tamanho
                somatorio += matriz[i,j] * vetor_sol[j]
            end

        vetor_sol[i] = (vetor[i]-somatorio)/matriz[i,i]

    end

    return vetor_sol

end

# -------------- FIM DA FUNÇÃO --------------

# ------------ ALGORTIMO DA DECOMPOSIÇÃO LU COM USO DO TEOREMA DA DECOMPOSIÇÃO LU.
"""
lu(A). ONDE "A" É UMA MATRIZ DE ORDEM n x n.

ESSA FUNÇÃO VAI DECOMPOR A MATRIZ "A" NO PRODUTO LU ONDE "L" É UMA MATRIZ TRIANGULAR INFERIOR E "U" É UMA MATRIZ TRIANGULAR SUPERIOR.

EX: A = [2 3 4; 2 -1 1; 2 2 6]

lu(A) ==> L = [1 0 0; 1 1 0; 1 0.25 1] e U = [2 3 4; 0 -4 -3; 0 0 2.75]

PARA TIRAR A PROVA REAL:

julia) L, U = lu(A)

julia) L*U == A

OU

julia) isapprox(L*U, A)

O RESULTADO FINAL SERÁ UM VALOR BOOLEANO "true" OU "false".
"""
function lu(matriz) # função responsável por verificar se é possivel decompor a matriz de coeficientes A em duas matrizes LU. Se possível, vai decompor a matriz A e mostrar as duas matrizes L e U obtidas.

    # criamos as variáveis tamanho, menor_principal e continuar
    # tamanho: recebe o tamanho da matriz quadrada inserida.
    # menor_principal: cria uma matriz cheia de zeros de tamanho igual a tamanho-1

    # ---------- VERIFICAÇÃO SE A MATRIZ INSERIDA É QUADRADA OU NÃO ----------------

    linha = size(matriz, 1)
    coluna = size(matriz, 2)

    if linha != coluna
        error("A matriz inserida não é quadrada e sim uma $linha x $coluna.")
    end

    # ---------- FIM DA VERIFICAÇÃO --------------

    tamanho = size(matriz, 1)
    matriz_l = Matrix{Float64}(I, tamanho, tamanho)
    matriz_u = zeros(Float64, tamanho, tamanho)
    println("Iniciando Decomposição...\n")

    # ------------- ALGORITMO DE DECOMPOSIÇÃO LU -------------.

    for i in 1:tamanho # iniciamos o processo como se estivessemos calculando um produto matricial usualmente, então escolhemos uma linha e realizamos o produto com todas as outras colunas.
        
        for j in 1:tamanho 

            if i <= j && i == 1 # o caso em que i <= j e i = 1 que é quando os coeficientes u_1j é igual a a_1j.

                matriz_u[i, j] = matriz[i, j]

                if matriz_u[i,i] == 0 
                    error("o elemento u_$i$i é nulo. Pelo Teorema da Decomposição LU, essa matriz não admite decomposição.")
                end

            elseif i <= j # segue para os próximos casos em que dependemos do valor l_ij. O algoritmo é igual a fórmula matemática feita.

                somatorio_u = 0.0 # definimos a somatorio_u aqui para resetar toda vez que começar uma nova conta na próxima linha.

                for k in 1:i-1

                    somatorio_u += matriz_l[i,k]*matriz_u[k,j]

                end

                matriz_u[i, j] = matriz[i, j] - somatorio_u

                if matriz_u[i,i] == 0
                    error("o elemento u_$i$i é nulo. Pelo Teorema da Decomposição LU, essa matriz não admite decomposição.")
                end

            elseif i > j # o caso em que a linha é maior que a coluna. Neste caso é calculado o valor de l_ij pois é possível isolá-lo. O algoritmo é literalmente igual a fórmula feita matemáticamente.

                somatorio_l = 0.0

                for k in 1:j-1

                    somatorio_l += matriz_l[i,k]*matriz_u[k,j]

                end

                matriz_l[i,j] = (matriz[i, j] - somatorio_l)/matriz_u[j, j]

            end
        end
    end

    

    # -------------- FIM DO ALGORITMO ------------------.

    # O programa finaliza e retorna as matrizes L e U cujo produto da a matriz A.

    produto_lu = matriz_l*matriz_u # salvamos o produto LU.

    if isapprox(norm(matriz - (produto_lu)), 0) == false # será se é uma boa comparação?

        erro_rel = norm(matriz - (produto_lu))/norm(matriz) # a norma padrão da julia é a norm p = 2 que é igual a norma encontrada em bibliografias de álgebra linear. A norma de Frobenius.
        println("O resultado do produto LU foi aproximado por um erro relativo de: $erro_rel\n")

    end

    println("---------Decomposição finalizada!---------\n")

    return matriz_l, matriz_u # retorno das matrizes L e U para poder trabalhar com a solução via LU do sistema linear desejado.

end

# --------------- FIM DA FUNÇÃO --------------------

# --------------- FUNÇÃO: Resolver o sistema linear Ax = b no formato (LU)x = b ------------------.

# É necessário inserir informações como: A matriz triangular inferior L, A matriz triangular superior U e o vetor de coeficientes b.
"""
resolver_lu(L, U ,b). A FUNÇÃO RESOLVE O SISTEMA LINEAR Ax = b PELO MÉTODO (LU)x = b.

"L" E "U" SÃO AS MATRIZES OBTIDAS PELA DECOMPOSIÇÃO LU E "b" É O VETOR DOS TERMOS INDEPENDENTES.

A FUNÇÃO UTILIZA AS FUNÇÕES sti() E sts() PARA O CÁLCULO.

"""
function resolver_lu(L, U, b)

    # -------------- VERIFICAÇÕES ANTES DE INICIAR O PROCESSO DE SOLUÇÃO DE SISTEMA POR LU --------------

    linha_b = size(b, 1)
    linha_L = size(L, 1)


    if linha_b != linha_L
        error("o vetor b inserido tem dimensão $linha_b enquanto que a matriz L tem dimensão $linha_L\n Portanto, não é possível realizar o produto L^(-1) * b")
    end


    # ATENÇÃO !
    # como estamos usando o determinante, acredito que seja mais eficiente usar apenas um checador se os elementos da diagonal principal não são nulos!

    if det(U) == 0 
        error("Como o determinante da matriz U é nula, não podemos realizar U^(-1) * y.\n Porntanto, impossível de resolver.")
    end

    # ---------------  FIM DAS VERIFICAÇÕES ----------------------
    # 
    # --------------- ALGORITMO PARA RESOLVER SISTEMA LINEAR VIA LU. UTILIZAMOS AS FUNÇÕES STI() E STS() JÁ CRIADAS PARA APROVEITAR ------------------
    
    y = sti(LowerTriangular(L), b) # obtemos o vetor y para usar depois na função sts()
    vetor_sol = sts(UpperTriangular(U), y) # aqui obtemos o vetor solução do sistema!

    println("O vetor solução é: $vetor_sol")

    return vetor_sol

end

# ---------------- FUNÇÃO: DECOMPOR A MATRIZ EM UM PRODUTO R^T * R (MÉTODO DE CHOLESKY) -----------------
"""
chol(A). ONDE A É UMA MATRIZ SIMÉTRICA DEFINIDA POSITIVA DE ORDEM n x n.
A FUNÇÃO VAI DECOMPOR A MATRIZ A NO PRODUTO R^T * R VIA MÉTODO DE CHOLESKY.
A MATRIZ "R^T" É A MATRIZ TRANSPOSTA DE "R" E "R" É UMA MATRIZ TRIANGULAR SUPERIOR.

EX: A = [5 2; 2 1]

chol(A) ==> R = [2.23607 0.89443; 0 0.44721] e R^T = [2.23607 0; 0.89443 0.44721]

PARA TIRAR A PROVA REAL:

julia) R, T = chol(A)
julia) T * R == A

OU

julia) isapprox(T*R, A)

O VALOR FINAL SERÁ UM VALOR BOOLEANO "true" OU "false".
"""
function chol(matriz)

    println("Iniciando Decomposição via Cholesky....\n")

    # PRIMEIRO CHECAMOS SE A MATRIZ INSERIDA É SIMÉTRICA OU NÃO. BASTA VERIFICAR SE A = A^T

    matriz_t = transpose(matriz)

    if isapprox(matriz_t, matriz) == false

        error("A matriz inserida não é simétrica! A ≠ A^T.")

    end

    # APÓS A VERIFICAÇÃO, INICIAMOS O CÁLCULO DA MATRIZ R VIA MÉTODO DE CHOLESKY.

    
    tamanho = size(matriz, 1) # PEGAMOS A DIMENSÃO USANDO COMO REFERÊNCIA O N° DE LINHAS DA MATRIZ.
    matriz_r = zeros(Float64, tamanho, tamanho)

    # INICIAMOS O ALGORITMO.

    for i in 1:tamanho 

        for j in i:tamanho

            if j == i && i == 1 # CASO i = j E i = 1, CASO EM QUE ESTAMOS CALCULANDO O ELEMENTO r_11.

                if matriz[i, i] < 0 || matriz[i, i] == 0 

                    r = matriz[i, i]

                    error("O termo r_$i$i = $r é um valor não positivo. Portanto, não é possível fazer a decomposição via cholesky.")
                
                else

                    matriz_r[i,i] = sqrt(matriz[i, i])

                end

            elseif j == i # CASO i = j ONDE O SOMATÓRIO APARECE. CASO EM QUE CALCULAMOS r_ii.

                somatorio_1 = 0.0

                for k in 1:i-1

                    somatorio_1 += (matriz_r[k, i])^2

                end

                result = matriz[i, i] - somatorio_1

                if result < 0 || result == 0

                    error("O termo r_$i$i = $result é não positivo. Portanto, não é possível fazer a decomposição via cholesky.")

                else

                    matriz_r[i, i] = sqrt(result)

                end

            else    # CASO EM QUE CALCULAMOS r_ij

                somatorio_2 = 0.0

                for k in 1:i-1

                    somatorio_2 += matriz_r[k, i]*matriz_r[k, j]

                end

                matriz_r[i, j] = (matriz[i, j] - somatorio_2)/matriz_r[i, i]

            end

        end

    end

    matriz_r_t = transpose(matriz_r)

    println("-------------- Decomposição via Cholesky finalizada! ---------------\n")

    return matriz_r, matriz_r_t

end

# ------------------ FIM DA FUNÇÃO ----------------------

# ----------------- FUNÇÃO: RESOLVER SISTEMA LINEAR VIA CHOLESKY. IMPLEMENTAR? --------------------

function resolver_chol(R, T, b)

end

# ---------------- FIM DA FUNÇÃO --------------------


# -------------------- FUNÇÃO QUE TESTA lu() para matrizes variando de dimensão 1 até 100. IMPLEMENTAR? ----------------------------

function testar_lu()
    
    println("Iniciando teste dos algortimos se estão funcionando corretamente!!! \n Quantas iterações deseja?")

    num_iteracoes = parse(Int64,readline())
    loop = 1

    while loop <= num_iteracoes

        println("$loop ° iteração.\n")
        # Vai ser escolhido um número de 1 até 100 para ser a dimensão da nossa matriz A.

        linha = rand(1:100)
        coluna = linha

        # Criaremos a matriz A com valores entre -100 a 100 escolhidos aleatoriamente.
        A = rand(-100.0:100.0, linha, coluna)

        # criaremos o vetor b com valores entre -100 a 100 escolhidos aleatoriamente.
        b = rand(-100.0:100.0, linha)

        # começaremos pela decomposição LU para obter as matrizes L e U.

        L, U = lu(A)

        # Após isso, faremos a solução via LU. caso não for possível ja vai ser retornado o erro da própria função.

        x = resolver_lu(L, U, b)

        erro = norm(b - ((L*U)*x))  # é calculado a diferença entre o vetor b e o vetor (LU)x que é justamente o b. Isso é feito em norma!

        print("o erro na $loop ° iteração foi de: $erro\n")
        loop += 1

    end

    println("O programa foi finalizado!\n")

end

# ---------------- FIM DA FUNÇÃO ---------------

# --------------- FUNÇÃO: MÉTODO ITERATIVO DE JACOBI-RICHARDSON ---------------------

# CRIAR UMA FUNÇÃO QUE VERIFICA SE A MATRIZ "A" INSERIDA É QUADRADA E SE O VETOR "b" É UM VETOR COLUNA?
"""
jacobi(A, b, x, tol). 

"A" É A MATRIZ QUADRADA DOS COEFICIENTES.
"b" É O VETOR DOS TERMOS INDEPENDENTES.
"x" É O VETOR SOLUÇÃO CHUTE INICIAL.
"tol" É A TOLERÂNCIA FIXA.

ESSA FUNÇÃO VAI CALCULAR UMA SOLUÇÃO APROXIMADA PARA O SISTEMA "Ax = b" BASEADO NOS ITENS INSERIDOS NA FUNÇÃO UTILIZANDO O MÉTODO ITERATIVO DE JACOBI.

EX: A = [2 1; 1 -2], b = [2, -2], x = [0, 0] e tol = 0.01

O RESULTADO VAI SER => [0.3984375, 1.1953125].
QUE É PRÓXIMO DA SOLUÇÃO EXATA DO SISTEMA ==> [0.4, 1.2]
"""
function jacobi(A, b, x, tol) 

    tamanho = size(A, 1) # OBTER O TAMANHO DA MATRIZ A
    x_iterativo = zeros(tamanho) # CRIAR UM VETOR "x aproximação k+1" PARA ARMAZENAR O RESULTADO DO MÉTODO ITERATIVO.
    H = zeros(tamanho, tamanho) # CONSTRUÇÃO DA MATRIZ ITERATIVA.
    g = zeros(tamanho) # CONSTRUÇÃODO VETOR DOS TERMOS INDEPENDENTES ITERATIVO.
    iteracao = 0 # PARA CONTABILIZAR A ITERAÇÃO
    loop = true # USAR NO LAÇO-WHILE.

    # VERIFICAR SE A MATRIZ "A" É DIAGONALMENTE DOMINANTE, USAREMOS A NORMA LINHA NOS ELEMETOS DA LINHA PARA COMPARAR COM O TERMO DA DIAGONAL PRINCIPAL DA MESMA LINHA.

    for i in 1:tamanho
    
        soma_1 = 0 # ARMAZENAR A SOMA DOS ELEMENTOS DA LINHA DESCONSIDERANDO O DA DIAGONAL PRINCIPAL.
    
        for j in 1:tamanho

            if j != i

                soma_1 += abs(A[i,j])

            end

            if abs(A[i,i]) < abs(soma_1) # VERIFICA SE O TERMO DA DIAGONAL PRINCIPAL É MENOR QUE A SOMA DOS OUTROS TERMOS DA MESMA LINHA.

                error("A matriz \"A\" não é diagonalmente dominante.")

            end
        end
    end

    println("A matriz \"A\" é diagonalmente dominante! Iniciando o Método de Jacobi\n")

    # APÓS A VERIFICAÇÃO, INICIAR O MÉTODO DE JACOBI. VAMOS CONSTRUIR A MATRIZ ITERATIVA E A MATRIZ ITERATIVA DOS TERMOS INDEPENDENTES.

    for i in 1:tamanho

        g[i] = b[i]/A[i, i] # CONSTRUÇÃO DO VETOR INDEPENDENTE ITERATIVO.

        for j in 1:tamanho

            if j != i # CONSTRUÇÃO DA MATRIZ ITERATIVA.
                
                H[i, j] = -A[i, j]/A[i, i]

            end

        end

    end

    # ENTRA NO LAÇO-WHILE E INICIA O ALGORITMO PARA CRIAR SOLUÇÕES QUE SE APROXIMAM DA SOLUÇÃO DO SISTEMA.

    while loop == true 

        iteracao += 1

        x_iterativo = H * x + g
        
        erro_relativo = norm(x_iterativo - x)/norm(x_iterativo) # CALCULA O ERRO RELATIVO PARA COMPARAR COM A TOLERANCIA FIXA.

        if erro_relativo < tol # SE O ERRO RELATIVO É MENOR QUE A TOLERANCIA, OBTEMOS A SOLUÇÃO APROXIMADA.

            loop = false

        else 

            x = x_iterativo # MUDANÇA DE "x" PARA RECEBER O RESULTADO E CONTINUAR A PRÓXIMA ITERAÇÃO.

        end

        

    end

    x_aproximado = x_iterativo

    println("Iteração concluida! n° de iterações: $iteracao")
    return x_aproximado

end

# ---------------- FIM DA FUNÇÃO -----------------------

# ---------------- FUNÇÃO: MÉTODO ITERATIVO DE GAUSS-SEIDEL ----------------------

function gauss()

end