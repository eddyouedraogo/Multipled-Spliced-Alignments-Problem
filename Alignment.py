''''Multipe Spliced Alignment Project'''
import os
'''La classe Alignments representera les alignments du segments alignees avec le gene'''
class Alignments:
    def __init__(self, idGene, idSeq, id, size, a, b, s=0, e=0, cons=0, exon="-"):
        self.id = id
        self.size = size
        self.a = a
        self.b = b
        self.s = s
        self.e = e
        self.cons = cons
        self.exon = exon
        self.idGene = idGene
        self.idSeq = idSeq


'''Le segment represente la sequence alignee contre le gene'''
class Segment:
    alignments = []  # Represents the aligned segments with the Genes 
    compatibles = []

    def __init__(self, id, size, cons, no):
        self.id = id
        self.size = size
        self.cons = cons
        self.no = no


'''La classe gene representes le Gene'''
class Genes:
    segments = []  # Represents the segments in the file
    compatibility = []

    def __init__(self, id):
        self.id = id

'La Classe Node represente les noeuds. ' \
'Lorsque que le gene a une portion debutant par s et finissant par e ' \
'a des sequences qui sont alignee a cette emplacement on creer un noeud et on rajoute ces sequences et leur informations'
class Node:
    children = []

    def __init__(self, idGene, s, e):
        self.idGene = idGene
        self.s = s
        self.e = e
'Apres avoir creer les noeuds on creer des blocks dans lesquels on mettra les genes qui ont les memes sequences ensembles.'
class Block:
    multiblocks = []
    
    def __init__(self,id):
        self.id = id
        
def readFile(filename):
    listOfGenes = []  # We Put all the genes that will be read into this list 
    geneAddedPos = 0  # To we can keep up with the genes position inside of the list
    alignmentList = []  # Hold the lines to temporary for us to extract the info and then use it later when we have all the info to create a segment 
    wordList = []  # List of words inside of the file 
    wordCount = 0  # Number of Words in the file
    emptycount = 0  # Helps count the empty lines
    lineCount = 0  # Count the number of lines
    linestate = ""  # Line state Empty, First or just Line 
    seq = Segment
    idAlignments = 0  # Alignments id
    gene = Genes
    align = Alignments
    posSegment = 0
    lC = 0
    listOfNodes = []
    file = open(filename, 'r')
    
    for line in file:

        if not line.split():
            linestate = "Empty Line"
            lineCount = 0  # lineCount nous permet de cibler la premiere ligne de la table 
            # Lorsqu'il y a une ligne vide on remet line count a zero
            # Ainsi lorsqu'on lit la premiere ligne on peut assigner les valeurs de Segments et Genes
        elif line.split():
            if lineCount == 0:
                linestate = "First Line"
            else:
                linestate = "Line"
            lineCount += 1

            for word in line.split():
                wordList.append(word)
                wordCount += 1

        if linestate == "First Line":
            lC += 1
            tempListSegment = wordList
            '''If the list is empty then we add elements otherwise we check if the genes is already here and then add it'''
            if len(listOfGenes) == 0 or len(listOfGenes) != 0 and wordList[1] != listOfGenes[len(listOfGenes) - 1].id:
                gene = Genes(wordList[1])  # The second word is always the id of the Gene
                listOfGenes.append(gene)

            geneAddedPos = len(listOfGenes) - 1  # Position of the gene in the list
            posSegment = len(gene.segments) - 1

        elif linestate == "Line":
            lC += 1
            if len(wordList) == 9:
                align = Alignments(wordList[1], wordList[0], idAlignments, wordList[2], wordList[3], wordList[4], wordList[5],
                                   wordList[6],
                                   wordList[7], wordList[8])
            else:
                align = Alignments(wordList[1], wordList[0], idAlignments, wordList[2], wordList[3], wordList[4])

            alignmentList.append(align)
            idAlignments += 1

        wordList = []
        
    
    for i in listOfGenes:
        emptyLst  = []
        for j in alignmentList:
            if i.id == j.idGene:
                emptyLst.append(j)
                i.segments = emptyLst
    
    node = None
    
    for i in listOfGenes:
        i.segments.sort(key=lambda x:(int(x.s)))
        for j in range(len(i.segments)-1):
            currentSeq = i.segments[j]
            nextSeq = i.segments[j+1]
            if int(currentSeq.s) - 50 <= int(nextSeq.s) <= int(currentSeq.s) + 50:
                if int(currentSeq.e) == int(nextSeq.e):
                    #Check if already in the list before adding it return Node to append children to
                    node = nodeExist(listOfNodes,Node(currentSeq.idGene,currentSeq.s,currentSeq.e))
                    node.children.append(currentSeq)
            elif int(currentSeq.e) - 50 <= int(nextSeq.e) <= int(currentSeq.e) + 50:
                if int(currentSeq.s) == int(nextSeq.s):
                    node = nodeExist(listOfNodes, Node(currentSeq.idGene, currentSeq.s, currentSeq.e))
                    node = nodeExist(listOfNodes, Node(currentSeq.idGene, currentSeq.s, currentSeq.e))
                    node.children.append(currentSeq)
            else: 
                node = nodeExist(listOfNodes,Node(currentSeq.idGene,currentSeq.s,currentSeq.e))

    
    for j in listOfNodes:
        tempChildren = []
        for i in range(len(j.children)):
            if j.children[i].idGene == j.idGene and j.children[i].s == j.s:
                tempChildren.append(j.children[i])
        j.children = tempChildren



    listOfNodes.sort(key=lambda x :(int(x.s)))


    blockList = []
    allBlocks = []
    blockID = 0
    for i in range(len(listOfNodes)-1):
        ab1 = int(listOfNodes[i].e) - int(listOfNodes[i].s)
        ab2 = int(listOfNodes[i+1].e) - int(listOfNodes[i+1].s)
        if(ab1- 15 <= ab2 <= ab1 + 15):
            blockExist(blockList,listOfNodes[i])
        else:
            block = Block(blockID)
            block.multiblocks = blockList
            allBlocks.append(block)
            blockList = []
            blockID += 1
            

    for i in allBlocks:
        sizeOfMultiblocks = len(i.multiblocks)
        for j in range(sizeOfMultiblocks):
            print(i.multiblocks[j].idGene, i.multiblocks[j].s, i.multiblocks[j].e)
        if sizeOfMultiblocks >= 1:
            for k in i.multiblocks[0].children:
                print(k.idSeq, k.a, k.b)

        print('\n')

    loneNodes = []
    for i in allBlocks:
        if len(i.multiblocks) == 1:
            loneNodes.append(i.multiblocks)

    return listOfGenes


def nodeExist(listOfNodes,node):
    found = False
    for n in listOfNodes:
        if n.idGene == node.idGene and n.s == node.s and n.e == node.e:
            found = True
            return n
    if not found:
        listOfNodes.append(node)
        return node
    
def blockExist(blockList,block):
    found = False
    for b in blockList:
        if b.idGene == block.idGene and b.s == block.s and b.e == block.e:
            found = True
            break
    if not found:
        blockList.append(block)

def compatibilityWithLoneNode(loneNode, nodes):
    compatible = False
    for i in range(len(loneNode.children)-2):
        for j in range(len(nodes.children)-2):
            if(len(nodes)>1):
                if(loneNode.children[i] == nodes.children[i]):
                    compatible = True
                    return compatible
    return compatible



filename = input("Please enter your filename : ")
print('\n')

if os.path.exists(filename):
   with open(filename, 'rb') as f:
       try:
           readFile(filename)
       except IOError as e:
           print(e)


