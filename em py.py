from numpy import array # para poder usar funcionalidades de matriz, vetor, etc.

# A pricípio, uma array é uma coleção n dimensional de objetos.
# n dimensional pois podemos vizualiza-lo como um vetor. isso é:
# se é 1-dimensional, é do R¹. Se é 2-dimensional é do R². Se é R³ é 3-dimensional.
# Acredito que seja possível fazer uma Array n por m dimensional que fazemos uma array matriz.

# matriz = np.array([[1, 2, 3], [4, 5, 6]])

# LEMBRAR QUE O PYTHON FAZ CONTAGEM COMEÇANDO EM 0 E INDO ATÉ N-1. ONDE N É A QUANTIDADE DE ELEMENTOS.

# print(matriz.shape) # o resultado vai ser (2, 3). Onde (i, j) i é o número de linhas e j é o numero de colunas.

# print(matriz[0][1]) # quero obter o dado da primeira linha na segunda coluna, então é o valor 2

# FAREMOS AGORA O PRODUTO INTERNO OU PRODUTO ESCALAR ENTRE DOIS VETORES DO R^n.


################################ PARTE PRÁTICA DO PROGRAMA #################################


# ---------- FUNÇÃO PRODUTO INTERNO OU PRODUTO ESCALAR DEFINIDA DE ACORDO COM O LIVRO. -------------
# FUNÇÃO CRIADA DE ACORDO COM O LIVRO MATRIX COMPUTATIONS - GOLUB.

def dot(x, y):
    print("\n")
    c = 0 # criar uma variável número que armazena esse resultado.

    linhas_x = x.shape # para receber o n° de colunas de a
    linhas_y = y.shape # para receber o n° de colunas de b.

    
    if linhas_x[0] != linhas_y[0]:

        print("não é possível realizar o produto! possuem dimensão diferente.\n")
        return 0
    
    else:

        print("Multiplicando os vetores:\n") # tem que colocar [0] pois o resultado retornado por shape sempre é uma tupla.
        print(f"{x} * {y}\n")

        for j in range(0, linhas_x[0]): # vamos desde a posição inicial até a final fazendo o produto e somando o resultado.

            c = c + x[j]*y[j] # pegamos os termos respectivos para fazer a multiplicação.

    print("Produto escalar realizado com sucesso!")
    print(f"Resultado: {c}\n")

    return c

# ---------------------------------------------------

# FUNÇÃO "Scalar a x plus y" (Saxpy). Nesse caso vamos fazer y[i] = y[i] + a*x[i]
# OBS: IMAGINO QUE É FEITO UMA SUBSTITUIÇÃO DO VETOR JA EXISTENTE PARA OCUPAR MENOS MEMÓRIA.

def saxpy(x, y, a):

    # Aqui obtemos os números de linhas e colunas dos vetores x e y.

    old_y = y # só para printar no resultado final a fim de verificar.

    linhas_x = x.shape
    linhas_y = y.shape

    

    if linhas_x[0] != linhas_y[0]:

        print("não é possivel realizar a soma! possuem dimensão diferente.\n")
        return 0
    
    else:

        for i in range(0, linhas_x[0]):

            y[i] = y[i] + a*x[i] # realizamos o produto a * x[i] e ai somamos com y[i].

    print(f"Soma realizada com sucesso entre {old_y} e {a*x}\n")
    print(f"Resultado:{y}\n")

    return y

# -----------------------------------------------------------------

# FUNÇÃO: GENERALIZED A times x plus y ou Gaxpy.
def gaxpy(x, y, A):

    linhas_x = x.shape # retorna o número de linhas do vetor x.
    linhas_y = y.shape # retorna o número de colunas do vetor y.

    dim_A = A.shape # para uma matriz, o shape retorna primeiro o número de linhas e depois o número de colunas.

    # baseado no algoritmo y[i] = y[i] + A[i, j]*x[i]
    # é preciso que i de x seja igual a j de A e que i de A seja igual a i de y.

    if linhas_x[0] != dim_A[1] or linhas_y[0] != dim_A[0]:

        print("Não é possível realizar o algoritmo gaxpy!\n")
        print("Ou o n° de linhas de x não bate com o n° de colunas de A\n")
        print("Ou o n° de linhas de y não bate com o n° de linhas de A\n")

        return 0

    else:

        for i in range(0, linhas_y[0]):

            for j in range(0, dim_A[1]):

                print(y[i], A[i, j], x[j])
                y[i] = y[i] + A[i, j]*x[j] # Só fazer no caderno para entender o raciocínio do que esta acontecendo aqui.

                


    print("Algoritmo Gaxpy realizado com sucesso!\n")
    print(f"Resultado: {y}")

    return y


vetor_x = array([2, 4])
vetor_y = array([1, 3])
matriz_A = array([[1, 1],[1, 1]])

# NÃO ESQUECER QUE ESTAMOS SUBSTITUINDO vetor_y POIS NÃO CRIAMOS OUTRO PARA OCUPAR MAIS MEMÓRIA.

dot(vetor_x, vetor_y) # retorna 100 ao fazer produto escalar de vetor_a com vetor_b

saxpy(vetor_x, vetor_y, 1) # retorna a soma entre vetor_x e vetor_y.

gaxpy(vetor_x, vetor_y, matriz_A) # RESULTADO FINAL É [7, 9]

    
