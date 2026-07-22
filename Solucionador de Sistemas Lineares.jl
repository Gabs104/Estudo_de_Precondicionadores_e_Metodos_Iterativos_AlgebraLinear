# -------- IMPORTAÇÕES DE BIBLIOTECAS --------

# A FUNÇÃO lu() FOI FEITA CONSIDERANDO RESTRIÇÃO DE QUE NÃO DEVE EXISTIR ZEROS NA DIAGONAL PRINCIPAL DA MATRIZ U DEVIDO O FATO DO OBJETIVO DO PROJETO DE BUSCAR RESOLVER SISTEMAS LINEARES.
# A FUNÇÃO jacobi() UTILIZA A VERIFICAÇÃO VIA ESTRITAMENTE DIAGONALMENTE DOMINANTE (abs(a_ii) > sum(a_ij). j in 1:tamanho e j ≠ i).
# A FUNÇÃO gauss() UTILIZA O CRITÉRIO DE SASSENFELD.
# A FUNÇÃO chol() EXIGE QUE A MATRIZ SEJA SIMÉTRICA E DEFINIDA POSITIVA (A = A^T e os elementos da diagonal principal são positivos).

# SUJEITO A MUDANÇAS.

using Random # IMPORTAR PARA USAR A BIBLIOTECA DE GERAÇÃO ALEATÓRIA.

# -------- FIM DAS IMPORTAÇÕES DE BIBLIOTECAS ----------
"""

FORNECE AS FUNÇÕES EXISTENTES NO PROGRAMA.

"""
function help()
    println("\nfunções existentes: sts(), sti(), lu(), chol(), jacobi(), gauss()\ndigite \"?\" e a função sem/com o parênteses para obter mais informações.")
end

function ismatrizquadrada(A) # FUNÇÃO QUE VERIFICA SE MATRIZ É QUADRADA.

    linha = size(A, 1)
    coluna = size(A, 2)

    if linha != coluna
        error("\nA matriz inserida não é quadrada e sim uma $linha x $coluna.")
    end
    
end

function isvetorcoluna(b) # FUNÇÃO QUE VERIFICA SE O VETOR DOS TERMOS INDEPENDENTES É COLUNA.

    if LinearAlgebra.ndims(b) != 1
        error("\no vetor dos termos independentes não é um vetor coluna!")
    end

end

function criterios(A, b, x, digitos, limite) # VERIFICA MAIS CONDIÇÕES PARA OS MÉTODOS ITERATIVOS.

    ismatrizquadrada(A)

    isvetorcoluna(b)

    isvetorcoluna(x)
    
    if digitos < 1 || typeof(digitos) == Char || typeof(digitos) == String || typeof(digitos) == Bool # VERIFICA SE É UM NÚMERO >= 1 OU NÃO.
        error("\no valor de n tem que ser um número maior ou igual a 1.")
    end

    if limite < 0
        error("o limite deve ser um número positivo ou nulo.")
    end

end



# ------------- ALGORITMO PARA RESOLVER SISTEMAS TRIANGULARES INFERIORES -------------------
"""

`sti(A, b).` 

`"A:Matrix"` É UMA MATRIZ TRIANGULAR INFERIOR DOS COEFICIENTES DE ORDEM n x n.

`"b::Vector"` É O VETOR DOS TERMOS INDEPENDENTES DE ORDEM n x 1.

ESSA FUNÇÃO VAI RESOVER UM SISTEMA LINEAR Ax = b USANDO UM SISTEMA TRIANGULAR INFERIOR.

# EX: ```A = [1 0 0; 2 3 0; 1 4 5] e b = [1, 3, 6]```

```sti(A, b) = [1, 0.3333333..., 0.73333333...]```

# PARA TIRAR A PROVA REAL:

```julia) x = sti(A, b)```

```julia) A*x```

```O RESULTADO SERÁ O VETOR b.```

"""
function sti(A, b) # o algoritmo de fato que resolve o sistema triangular inferior pelo algoritmo descoberto no livro: "Cálculo Numérico: aprendizagem com Apoio de Software"

    ismatrizquadrada(A)  # VERIFICA SE A MATRIZ É QUADRADA.

    isvetorcoluna(b) # VERIFICA SE O VETOR É COLUNA.

    # Toda vez que iniciarmos essa função definimos o vetor solução e o tamanho da matriz.

    tamanho = size(A, 1) # fornece a o n° de linhas da matriz, como é quadrada o numero de colunas é igual.

    vetor_sol = zeros(Float64, tamanho) # cria um vetor cheio de zeros com tamanho igual ao número de linhas da matriz que se esta resolvendo. Esse é o vetor que vai armazenar a solução.

    for i in 1:tamanho # fazemos um laço for que depende do tamanho da matriz.

        somatorio = 0.0 # para cada iteração, a soma linear muda.

        for j in 1:i-1 # é o somatório de j indo de 1 até i-1 somando a_ij*x_j
            somatorio += A[i, j] * vetor_sol[j]
        end

        # aqui calculamos o valor de x após ter realizado o somatório e armazenamos no vetor solução.
        vetor_sol[i] = (b[i] - somatorio) / A[i, i]

    end

    return vetor_sol

end

# --------------- FIM DA FUNÇÃO ---------------


# ------------------- ALGORITMO PARA RESOLVER SISTEMAS TRIANGULARES SUPERIORES -------------------------

"""

`sts(A, b).`

`"A::Matrix"` É UMA MATRIZ TRIANGULAR SUPERIOR DOS COEFICIENTES DE ORDEM n x n. 

`"b::Vector"` É O VETOR DOS TERMOS INDEPENDENTES.

ESSA FUNÇÃO VAI RESOLVER O SISTEMA LINEAR `Ax = b` VIA SISTEMA TRIANGULAR SUPERIOR.

# EX: ```A = [1 3 4; 0 4 3; 0 0 -1] e b = [3, -2, 1]```

```sts(A, b) = [6.25, 0.25, -1.0]```

```PARA TIRAR A PROVA REAL:```

```julia) x = sts(A, b)```

```julia) A*x```

# O RESULTADO SERÁ O VETOR "b".

"""
function sts(A, b) # o algoritmo de fato que resolve o sistema triangular superior pelo algoritmo descoberto no livro: "Cálculo Numérico: aprendizagem com Apoio de Software"

    ismatrizquadrada(A)  # VERIFICA SE A MATRIZ É QUADRADA.

    isvetorcoluna(b) # VERIFICA SE O VETOR É COLUNA.

    # toda vez que iniciamos essa função definimos a variável "tamanho" que recebe o tamanho da matriz e a variável "vetor_sol" que recebe o nosso vetor x solução.

    tamanho = size(A, 1) # o size retorna quantas linhas tem essa matriz, se colocar o 2 ele retorna a quantidade de colunas.

    vetor_sol = zeros(Float64, tamanho) # ele é composto com entradas de zeros que são do tipo float64.

    for i in tamanho:-1:1 # laço for que vai diminuindo desde a n-ésima linha até a primeira linha.

        somatorio = 0.0 # a soma linear que

            for j in i+1:tamanho
                somatorio += A[i,j] * vetor_sol[j]
            end

        vetor_sol[i] = (b[i]-somatorio)/A[i,i]

    end

    return vetor_sol

end

# -------------- FIM DA FUNÇÃO --------------

# ------------ ALGORTIMO DA DECOMPOSIÇÃO LU COM USO DO TEOREMA DA DECOMPOSIÇÃO LU.
"""

lu(A). 

`"A::Matrix"` É UMA MATRIZ DE ORDEM n x n.

ESSA FUNÇÃO VAI DECOMPOR A MATRIZ "A" NO PRODUTO LU ONDE "L" É UMA MATRIZ TRIANGULAR INFERIOR E "U" É UMA MATRIZ TRIANGULAR SUPERIOR.

# EX: A = [2 3 4; 2 -1 1; 2 2 6]

lu(A) ==> `L = [1 0 0; 1 1 0; 1 0.25 1]` e `U = [2 3 4; 0 -4 -3; 0 0 2.75]`

"""
function lu(A) # função responsável por verificar se é possivel decompor a matriz de coeficientes A em duas matrizes LU. Se possível, vai decompor a matriz A e mostrar as duas matrizes L e U obtidas.

    # criamos as variáveis tamanho, menor_principal e continuar
    # tamanho: recebe o tamanho da matriz quadrada inserida.
    # menor_principal: cria uma matriz cheia de zeros de tamanho igual a tamanho-1

    ismatrizquadrada(A) # VERIFICA SE MATRIZ É QUADRADA OU NÃO.

    tamanho = size(A, 1)

    L = Matrix{Float64}(I, tamanho, tamanho)
    U = zeros(Float64, tamanho, tamanho)

    println("Iniciando Decomposição...\n")

    # ------------- ALGORITMO DE DECOMPOSIÇÃO LU -------------.

    for i in 1:tamanho # iniciamos o processo como se estivessemos calculando um produto matricial usualmente, então escolhemos uma linha e realizamos o produto com todas as outras colunas.
        
        for j in 1:tamanho 

            if i <= j # segue para os próximos casos em que dependemos do valor l_ij. O algoritmo é igual a fórmula matemática feita.

                somatorio_u = 0.0 # definimos a somatorio_u aqui para resetar toda vez que começar uma nova conta na próxima linha.

                for k in 1:i-1 # 

                    somatorio_u += L[i,k]*U[k,j]

                end

                U[i, j] = A[i, j] - somatorio_u 

                if U[i,i] == 0
                    error("o elemento u_$i$i é nulo. Pelo Teorema da Decomposição LU, essa matriz não admite decomposição.")
                end

            elseif i > j # o caso em que a linha é maior que a coluna. Neste caso é calculado o valor de l_ij pois é possível isolá-lo. O algoritmo é literalmente igual a fórmula feita matemáticamente.

                somatorio_l = 0.0

                for k in 1:j-1 

                    somatorio_l += L[i,k]*U[k,j]

                end

                L[i,j] = (A[i, j] - somatorio_l)/U[j, j] 

            end
        end
    end

    

    # -------------- FIM DO ALGORITMO ------------------.

    # O programa finaliza e retorna as matrizes L e U cujo produto da a matriz A.

    produto_lu = L*U # salvamos o produto LU.

    if isapprox(LinearAlgebra.norm(A - (produto_lu)), 0) == false 

        erro_rel = LinearAlgebra.norm(A - (produto_lu))/LinearAlgebra.norm(A) # a norma padrão da julia é a norm p = 2 que é igual a norma encontrada em bibliografias de álgebra linear. A norma de Frobenius.
        println("O resultado do produto LU foi aproximado por um erro relativo de: $erro_rel\n")

    end

    println("---------Decomposição finalizada!---------\n")

    return L, U # retorno das matrizes L e U para poder trabalhar com a solução via LU do sistema linear desejado.

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


    isvetorcoluna(b) # VERIFICA SE O VETOR b É UM VETOR COLUNA.

    if linha_b != linha_L # vetificar se as linhas são iguais 
        error("o vetor b inserido tem dimensão $linha_b enquanto que a matriz L tem dimensão $linha_L\n Portanto, não é possível realizar o produto L^(-1) * b")
    end


    # ATENÇÃO !
    # como estamos usando o determinante, acredito que seja mais eficiente usar apenas um checador se os elementos da diagonal principal não são nulos!

    if LinearAlgebra.det(U) == 0 
        error("Como o determinante da matriz U é nula, não podemos realizar U^(-1) * y.\n Porntanto, impossível de resolver.")
    end

    # ---------------  FIM DAS VERIFICAÇÕES ----------------------
    # 
    # --------------- ALGORITMO PARA RESOLVER SISTEMA LINEAR VIA LU. UTILIZAMOS AS FUNÇÕES STI() E STS() JÁ CRIADAS PARA APROVEITAR ------------------
    
    y = sti(LinearAlgebra.LowerTriangular(L), b) # obtemos o vetor y para usar depois na função sts()
    vetor_sol = sts(LinearAlgebra.UpperTriangular(U), y) # aqui obtemos o vetor solução do sistema!

    println("O vetor solução é: $vetor_sol")

    return vetor_sol

end

# ---------------- FUNÇÃO: DECOMPOR A MATRIZ EM UM PRODUTO R^T * R (MÉTODO DE CHOLESKY) -----------------
"""

`chol(A).

`"A::Matrix"` É UMA MATRIZ SIMÉTRICA DEFINIDA POSITIVA DE ORDEM n x n.

A FUNÇÃO VAI DECOMPOR A MATRIZ A NO PRODUTO `"R^T * R"` VIA MÉTODO DE CHOLESKY.

A MATRIZ "R^T" É A MATRIZ TRANSPOSTA DE `"R"` E `"R"` É UMA MATRIZ TRIANGULAR SUPERIOR.

# EX: `A = [5 2; 2 1]`

# chol(A) ==> `R = [2.23607 0.89443; 0 0.44721] e R^T = [2.23607 0; 0.89443 0.44721]`

"""
function chol(A)

    ismatrizquadrada(A)

    println("Iniciando Decomposição via Cholesky....\n")

    # PRIMEIRO CHECAMOS SE A MATRIZ INSERIDA É SIMÉTRICA OU NÃO. BASTA VERIFICAR SE A = A^T

    if LinearAlgebra.issymmetric(A) == false

        error("A matriz inserida não é simétrica! A ≠ A^T.")

    end

    # APÓS A VERIFICAÇÃO, INICIAMOS O CÁLCULO DA MATRIZ R VIA MÉTODO DE CHOLESKY.

    
    tamanho = size(A, 1) # PEGAMOS A DIMENSÃO USANDO COMO REFERÊNCIA O N° DE LINHAS DA MATRIZ.
    R = zeros(Float64, tamanho, tamanho)

    # INICIAMOS O ALGORITMO.

    for i in 1:tamanho 

        for j in i:tamanho

            if j == i && i == 1 # CASO i = j E i = 1, CASO EM QUE ESTAMOS CALCULANDO O ELEMENTO r_11.

                if A[i, i] <= 0  

                    r = A[i, i]

                    error("O termo r_$i$i = $r é um valor não positivo. Portanto, não é possível fazer a decomposição via cholesky.")
                
                else

                    R[i,i] = sqrt(A[i, i])

                end

            elseif j == i # CASO EM QUE CALCULAMOS r_ii.

                somatorio_1 = 0.0

                for k in 1:i-1

                    somatorio_1 += (R[k, i])^2

                end

                result = A[i, i] - somatorio_1

                if result <= 0 

                    error("O termo r_$i$i = $result é não positivo. Portanto, não é possível fazer a decomposição via cholesky.")

                else

                    R[i, i] = sqrt(result)

                end

            else    # CASO EM QUE CALCULAMOS r_ij

                somatorio_2 = 0.0

                for k in 1:i-1

                    somatorio_2 += R[k, i]*R[k, j]

                end

                R[i, j] = (A[i, j] - somatorio_2)/R[i, i]

            end

        end

    end

    R_transposta = LinearAlgebra.transpose(R)

    println("-------------- Decomposição via Cholesky finalizada! ---------------\n")

    return R, R_transposta

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

        erro = LinearAlgebra.norm(b - ((L*U)*x))  # é calculado a diferença entre o vetor b e o vetor (LU)x que é justamente o b. Isso é feito em norma!

        print("o erro na $loop ° iteração foi de: $erro\n")
        loop += 1

    end

    println("O programa foi finalizado!\n")

end

# ---------------- FIM DA FUNÇÃO ---------------

# --------------- FUNÇÃO: MÉTODO ITERATIVO DE JACOBI-RICHARDSON ---------------------

# CRIAR UMA FUNÇÃO QUE VERIFICA SE A MATRIZ "A" INSERIDA É QUADRADA E SE O VETOR "b" É UM VETOR COLUNA?
"""

`jacobi(A, b, x, tol, digitos, lim).` 

# ESSA FUNÇÃO VAI CALCULAR UMA SOLUÇÃO APROXIMADA PARA O SISTEMA "Ax = b" BASEADO NOS ITENS INSERIDOS NA FUNÇÃO UTILIZANDO O MÉTODO ITERATIVO DE JACOBI.

`"A::Matrix"` É A MATRIZ QUADRADA DOS COEFICIENTES.

`"b::Vector"` É O VETOR DOS TERMOS INDEPENDENTES.

`"chute_inicial"` É O VETOR SOLUÇÃO CHUTE INICIAL.

`"tol"` É A TOLERÂNCIA FIXA.

`"digitos"` É QUANTIDADE DE CASAS DECIMAIS A SE CONSIDERAR DO RESULTADO FINAL DA ITERAÇÃO.

`"lim"` É O NÚMERO ITERAÇÕES A SER REALIZADA. CASO limite = 0 SERÁ FEITO ITERAÇÕES ATÉ PASSAR TOLERÂNCIA FIXA.

É FEITO A VERIFICAÇÃO DA MATRIZ SER DIAGONALMENTE DOMINANTE OU NÃO VIA A NORMA COLUNA.

# EX: ```A = [2 1; 1 -2], b = [2, -2], x = [0, 0], tol = 0.01 e digitos = 2, limite = 0```

O RESULTADO VAI SER => [0.4, 1.2].

"""
function jacobi(A, b, chute_inicial, tol, digitos, limite)

    x .= Vector{Float64}(chute_inicial)

    criterios(A, b, x, digitos, limite) # VERIFICAÇÃO DOS CRITÉRIOS PARA USAR O MÉTODO.

    tamanho = size(A, 1) # OBTER O TAMANHO DA MATRIZ A.

    x_iterativo = zeros(tamanho) # CRIAR UM VETOR "x aproximação k+1" PARA ARMAZENAR O RESULTADO DO MÉTODO ITERATIVO.

    H = zeros(tamanho, tamanho) # CONSTRUÇÃO DA MATRIZ ITERATIVA.

    g = zeros(tamanho) # CONSTRUÇÃODO VETOR DOS TERMOS INDEPENDENTES ITERATIVO.

    iteracao = 0 # PARA CONTABILIZAR A ITERAÇÃO

    loop = true # USAR NO LAÇO-WHILE.

    verificacao = true

    # VERIFICAR SE A MATRIZ "A" É DIAGONALMENTE DOMINANTE, USAREMOS A NORMA LINHA NOS ELEMETOS DA LINHA PARA COMPARAR COM O TERMO DA DIAGONAL PRINCIPAL DA MESMA LINHA.

    for i in 1:tamanho
    
        soma_1 = 0 # ARMAZENAR A SOMA DOS ELEMENTOS DA LINHA DESCONSIDERANDO O DA DIAGONAL PRINCIPAL.
    
        for j in 1:tamanho # 

            if j != i

                soma_1 += abs(A[i,j])

            end

            if abs(A[i,i]) <= abs(soma_1) && verificacao == true # VERIFICA SE O TERMO DA DIAGONAL PRINCIPAL É MENOR QUE A SOMA DOS OUTROS TERMOS DA MESMA LINHA. ESTRITAMENTE DIAGONALMENTE DOMINANTE.

                @warn "A matriz \"A\" não é diagonalmente dominante. Portanto, a convergência pode não ocorrer!"
                verificacao = false

                # SE APOS VERIFICA ISSO AQUI, SERIA BOM SE ESSA PARTE JA FOSSE PULADA PARA A ITERAÇÃO, CERTO?
            end
            

        end

    end

    if verificacao == true

        println("A matriz \"A\" é diagonalmente dominante! Iniciando o Método de Jacobi\n")

    end

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
        
        erro_relativo = LinearAlgebra.norm(x_iterativo - x)/LinearAlgebra.norm(x_iterativo) # CALCULA O ERRO RELATIVO PARA COMPARAR COM A TOLERANCIA FIXA.

        if any(isnan, x_iterativo) # VERIFICAÇÃO SE O VETOR ITERATIVO EXPLODIU PARA +∞ ou -∞

            error("A iteração levou $iteracao passos e divergiu da solução...\nTente precondicionar a matriz.")

        elseif erro_relativo < tol && limite == 0 # SE O ERRO RELATIVO É MENOR QUE A TOLERANCIA, OBTEMOS A SOLUÇÃO APROXIMADA.

            loop = false
            
        elseif iteracao == limite # CASO O USUÁRIO DESEJA APENAS FOCAR NO NÚMERO DE ITERAÇÕES.

            loop = false

        else
            x .= x_iterativo # MUDANÇA DE "x" PARA RECEBER O RESULTADO E CONTINUAR A PRÓXIMA ITERAÇÃO.
        end
             
    end

    
    x_aproximado = round.(x_iterativo, digits = digitos)
    println("Iteração concluida! n° de iterações: $iteracao")

    return x_aproximado

end

# ---------------- FIM DA FUNÇÃO -----------------------

# ---------------- FUNÇÃO: MÉTODO ITERATIVO DE GAUSS-SEIDEL ----------------------
"""

`gauss(A, b, chute_inicial, tol, digitos, limite).`

`- "A::Matrix"` É A MATRIZ DOS COEFICIENTES DO SISTEMA LINEAR.

`- "b::Vector"` É O VETOR DOS TERMOS INDEPENDENTES.

`- "chute_inicial::Vector"` É O VETOR CHUTE INICIAL.

`- "tol::Int64"` É A TOLERÂNCIA FIXA.

`- "digitos::Int64"` É QUANTIDADE DE CASAS DECIMAIS A SE CONSIDERAR DO RESULTADO FINAL DA ITERAÇÃO.

`- "limite::Int64"` É O NÚMERO LIMITE DE ITERAÇÕES A SER REALIZADA. 

CASO `- "limite = 0"` SERÁ FEITO ITERAÇÕES ATÉ PASSAR TOLERÂNCIA FIXA.

A FUNÇÃO VAI CALCULAR UMA APROXIMAÇÃO PARA A SOLUÇÃO DO SISTEMA LINEAR PELO
MÉTODO ITERATIVO DE GAUSS-SIEDEL. UTILIZA-SE O CRITÉRIO DE SASSENFELD.
É POSSÍVEL ARMAZENAR O VETOR SOLUÇÃO.

# EX: ```A = [10 2 1; 1 5 1; 2 3 10], b = [14, 11, 8], x = [0, 0, 0], ϵ = 0.01 e n = 1.``` 
 
`julia> gauss(A, b, x, tol)`

`julia> Iteração concluida em 4 passos!`

`julia> A solução do sistema é: [1, 2, 0]`

"""
function gauss(A, b, chute_inicial, tol, digitos, limite)

    x .= Vector{Float64}(chute_inicial) # APENAS PARA DIZER QUE O VETOR "x" ACEITA VALORES PONTOS FLUTUANTES.

    criterios(A, b, x, digitos, limite) # VERIFICAÇÃO DOS CRITÉRIOS PARA USAR O MÉTODO.

    tamanho = size(A, 1) # DIMENSÃO DA MATRIZ INSERIDA

    β = zeros(tamanho) # CRIAÇÃO DO VETOR BETA PARA VERIFICAÇÃO DO CRITÉRIO DE SASSENFIELD.

    x_iterativo = zeros(tamanho) # CRIAÇÃO DO VETOR PARA ITERAÇÃO (K+1).

    loop = true # PARA CONTINUAR O MÉTODO ITERATIVO.

    iteracao = 0 # PARA CONTABILIZAR AS ITERAÇÕES.

    # INICIAREMOS AQUI O CÁLCULO DO CRITÉRIO DE SASSENFIELD.

    for i in 1:tamanho

        soma_1 = 0 # PARA O PRIMEIRO SOMATÓRIO.
        soma_2 = 0 # PARA O SEGUNDO SOMATÓRIO.

        for j in 1:i-1 # PRIMEIRO SOMATÓRIO.

            soma_1 += abs(A[i, j]/A[i, i]) * β[j]

        end

        for j in i+1:tamanho # SEGUNDO SOMATÓRIO.

            soma_2 += abs(A[i, j]/A[i, i])

        end

        β[i] = soma_1 + soma_2 # INSERIMOS O RESULTADO NO VETOR BETA PARA ARMAZENAR.

    end

    if maximum(β) >= 1 # O CRITÉRIO EXIGE QUE O MÁXIMO DAS ENTRADAS DE BETA SEJA INFERIOR A 1.

        @warn "A matriz inserida não atende o critério de Sassenfeld!\n Portanto, a convergência pode não ocorrer!\n"

    else

        println("A matriz inserida atende o critério de sassenfeld! Portanto, a convergência ocorre rapidamente...\n")

    end

    # INICIAMOS AQUI O MÉTODO ITERATIVO

    while loop == true

        iteracao += 1

        for i in 1:tamanho

            soma_3 = 0 # PARA ARMAZENAR O RESULTADO DO PRIMEIRO SOMATÓRIO.
            soma_4 = 0 # PARA ARMAZENAR O RESULTADO DO SEGUNDO SOMATÓRIO.

            for j in 1:i-1 # PRIMEIRO SOMATÓRIO

                soma_3 += A[i, j]/A[i, i] * x_iterativo[j]

            end

            for j in i+1:tamanho # SEGUNDO SOMATÓRIO

                soma_4 += A[i, j]/A[i, i] * x[j]

            end

            x_iterativo[i] = b[i]/A[i, i] - soma_3 - soma_4 # ARMAZENAR O RESULTADO NO VETOR ITERAÇÃO.

        end
        
        erro_relativo = LinearAlgebra.norm(x_iterativo - x)/LinearAlgebra.norm(x_iterativo) # CALCULA O ERRO RELATIVO.

        # CASO A MATRIZ NÃO ATENDER O CRITÉRIO DE SASSENFELD, PODE OCORRER DE A ITERAÇÃO NÃO CONVERGIR PARA A SOLUÇÃO DO SISTEMA.
        # NESSE CASO, SERA ADICIONADO UM CONDICIONAL IF PARA CANCELAR A ITERAÇÃO E NÃO FICAR RODANDO O PROGRAMA INFINITAMENTE.
        # UMA MOTIVAÇÃO PARA O ESTUDO DE MÉTODOS DE PRECONDICIONAMENTO DE MATRIZES!
        # OBSERVAÇÃO: TESTE PARA MATRIZES DE DIMENSÃO MAIOR E VERIFIQUE QUE O CRITÉRIO NÃO É SATISFEITO FACILMENTE, NISSO O VETOR "x_iterativo" EXPLODE PARA +∞ ou -∞

        if any(isnan, x_iterativo)
            error("A iteração levou $iteracao passos e divergiu da solução...\nTente precondicionar a matriz.")
        end

        if erro_relativo < tol && limite == 0

            println(erro_relativo, tol)
            x_aproximado = round.(x_iterativo, digits = digitos)
            println("Iteração concluida em $iteracao passos!\n")
            println("A solução do sistema é: $x_aproximado\n")
            loop = false
            return x_aproximado

        elseif iteracao == limite

            x_aproximado = round.(x_iterativo, digits = digitos)
            println("Iteração concluida em $iteracao passos!\n")
            println("A solução do sistema é: $x_aproximado\n")
            loop = false
            return x_aproximado

        else

            x .= x_iterativo

        end

    end

end

# ------------- FIM DA FUNÇÃO ---------------

# ------------ FUNÇÃO: MÉTODO DOS GRADIENTES CONJUGADOS ---------------

"""
`gradconj(A, b, chute_inicial, tol, digitos).`

`- "A::Matrix"` É A MATRIZ DOS COEFICIENTES DO SISTEMA LINEAR.

`- "b::Vector"` É O VETOR DOS TERMOS INDEPENDENTES.

`- "chute_inicial::Vector"` É O VETOR CHUTE INICIAL.

`- "tol:Int64"` É A TOLERÂNCIA FIXA.

`- "digitos:Int64"` É QUANTIDADE DE CASAS DECIMAIS A SE CONSIDERAR DO RESULTADO FINAL DA ITERAÇÃO.

ESSA FUNÇÃO NÃO CONTÉM UM MÉTODO DE PRE-CONDICIONAMENTO EMBUTIDO.

"""
function gradconj(A, b, chute_inicial, tol, digitos)

    x .= chute_inicial  

    # VERIFICAÇÃO DAS ENTRADAS.

    ismatrizquadrada(A)
    isvetorcoluna(x)
    isvetorcoluna(b)

    # ANTES DE USAR O MÉTODO ITERATIVO DOS GRADIENTES CONJUGADOS, PRECISAMOS VERIFICAR SE>
    # "A" É UMA MATRIZ DEFINIDA POSITIVA.

    chol(A) # BASTA USARMOS A FUNÇÃO chol(A). SE ISSO FOR REALIZADO COM SUCESSO A MATRIZ "A" É SIMÉTRICA E DEFINIDA POSITIVA.

    tamanho = size(A, 1)

    iteracao = 0

    loop = true

    # PRIMEIRO CASO DA ITERAÇÃO ONDE CALCULAMOS "r" E DIZEMOS QUE "v" É IGUAL A "r".

    r = zeros(tamanho)
    v = zeros(tamanho)

    r .= b - A*x
    v .= r

    x_iterativo = zeros(tamanho)
    r_iterativo = zeros(tamanho)
    v_iterativo = zeros(tamanho)

    # SE A MATRIZ "A" INSERIDA É SIMÉTRIA E DEFINIDA POSITIVAMENTE, ENTÃO ELA ADMITE UMA DECOMPOSIÇÃO CHOLESKY.

    # APÓS ESSAS VERIFICAÇÕES, COMEÇAMOS O AQUI O MÉTODO DE GRADIENTES CONJUGADOS

    while loop == true
        
        iteracao += 1

        t = dot(r, r)/dot(v, A*v) # CÁLCULO DO COEFICIENTE "t".

        x_iterativo .= x + t*v # CÁLCULO DA PRÓXIMA APROXIMAÇÃO.

        r_iterativo .= r - t*(A*v) # CÁLCULO DO NOVO VETOR RESÍDUO.

        s = dot(r_iterativo,r_iterativo)/dot(r, r) # CÁLCULO DO COEFICIENTE "s".

        v_iterativo .= r_iterativo + s*v # CÁLCULO DO NOVO VETOR DIREÇÃO.

        x .= x_iterativo
        v .= v_iterativo

        if LinearAlgebra.norm(r) < tol || LinearAlgebra.norm(r) == 0 # VERIFICAMOS SE ESTÁ ABAIXO DA TOLERÂNCIA FIXA OU SE A NORMA É ZERO. SE SIM, OBTEMOS A SOLUÇÃO.

            println("Iteração concluida! n° de passos: $iteracao\n")

            x_solução = round.(x_iterativo, digits = digitos)

            loop = false

            return x_solução

        else

            r .= r_iterativo

        end

    end
end