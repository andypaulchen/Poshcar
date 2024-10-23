# Connect.py: Graph analytics on crystals

# 2 February 2024: I popped the question to Sophia last friday at Sungei Buloh. 
# The witnesses were 30 or so migratory birds at the tidal pools. Everyone was very happy that day

from poshcar.distance import * # Distances package
import networkx as nx # We are using networkx for graph analytics!!!!!
from pyvis.network import Network as net
from IPython.display import display, HTML

Image_colors = ['black', 'blue', 'blue', 'green', 'cyan', 'darkcyan', 'green', 'darkcyan', 'cyan', 'red', 'magenta', 'darkmagenta', 'yellow', 'peru', 'sienna', 'gold', 'dimgrey', 'rosybrown', 'red', 'darkmagenta', 'magenta', 'gold', 'rosybrown', 'dimgrey', 'yellow', 'sienna', 'peru']

#bar = u'\u0331' # bar below combining symbol
bar = '' # For undirected graph; reserve this variable in case needed
Image_labels = ['', 'c', 'c'+bar, 'b', 'bc', 'bc'+bar, 'b'+bar, 'b'+bar+'c', 'b'+bar+'c'+bar, 'a', 'ac', 'ac'+bar, 'ab', 'abc', 'abc'+bar, 'ab'+bar, 'ab'+bar+'c', 'ab'+bar+'c'+bar, 'a'+bar, 'a'+bar+'c', 'a'+bar+'c'+bar, 'a'+bar+'b', 'a'+bar+'bc', 'a'+bar+'bc'+bar, 'a'+bar+'b'+bar, 'a'+bar+'b'+bar+'c', 'a'+bar+'b'+bar+'c'+bar]

def isConnected(data, tolerance):
    # Input crystal structure and tolerance, returns sets of connected components
    # Collect headers
    atomspp = elemindices(data)
    res_indexed = atomspp['POSCAR Site']
    ns = len(res_indexed)
    bm = matrix_bonding(data, tolerance, verbose = False) # get bonding matrices
    
    G = nx.MultiGraph()
    for virtual in range(27):
        for i in range(ns):
            for j in range(ns)[i:]: 
                ethick = 1 if virtual == 0 else 5
                if bm[virtual][i][j] == 1: 
                    G.add_edge(res_indexed[i], res_indexed[j], weight=ethick, \
                               color=Image_colors[virtual], label=Image_labels[virtual], arrows = "")
                    
    connected_sets = sorted(nx.connected_components(G), key=len)
    return connected_sets, G

def Coordination_Graph(data, tolerance):
    # Return graph of coordinations, averaged over element
    runiq, bme = matrix_bonding_average(data, 'element', tolerance, verbose = True) # get bonding matrices
    
    G = nx.MultiDiGraph()
    for i in range(len(runiq)):
        for j in range(len(runiq)):
            if bme[i][j] > 0: 
                G.add_edge(runiq[i], runiq[j], weight=bme[i][j], label=str(bme[i][j]))
    
    return G

def Draw_Graph(G, output = "untitled_graph.html"):
    # Draws a graph
    netcryst = net(notebook=True, directed=True, cdn_resources = 'in_line')
    netcryst.from_nx(G)
    netcryst.set_options('''
    var options = {
    "edges": {
    "arrowStrikethrough": false,
    "color": {
    "inherit": false
    },
    "font": {
    "size": 10,
    "align": "top"
    },
    "smooth": true
    }
    }
    ''')

    html = netcryst.generate_html()
    with open(output, mode='w', encoding='utf-8') as fp:
        fp.write(html)
    display(HTML(html))