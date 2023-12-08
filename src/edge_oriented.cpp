#include "edge_oriented.h"
#include <set>
#include <algorithm>
#include <unordered_map>
#include "set_operation.h"

extern const int K, L;
extern unsigned long long N;

EBBkC_Graph_t::EBBkC_Graph_t() = default;

EBBkC_Graph_t::~EBBkC_Graph_t() {
    int i;

    delete [] edges;

    if (T) {
        for (i = 0; i < e_size; i++) delete [] T[i];
        delete [] T;
    }

    if (C) {
        for (i = 0; i < e_size; i++) delete [] C[i];
        delete [] C;
    }

    delete [] T_size;

    delete [] C_size;

    if (sub_v) {
        for (i = 0; i <= K; i++) delete [] sub_v[i];
        delete [] sub_v;
    }

    if (sub_e) {
        for (i = 0; i <= K; i++) delete [] sub_e[i];
        delete [] sub_e;
    }

    delete [] sub_v_size;

    delete [] sub_e_size;

    delete [] lab;

    if (DAG_deg) {
        for (i = 0; i <= K; i++) delete [] DAG_deg[i];
        delete [] DAG_deg;

    }

    if (G_deg) {
        for (i = 0; i <= K; i++) delete [] G_deg[i];
        delete [] G_deg;
    }


    delete [] col;

    if (DAG_adj) {
        for (i = 0; i < v_size; i++) delete [] DAG_adj[i];
        delete [] DAG_adj;
    }

    if (G_adj) {
        for (i = 0; i < v_size; i++) delete [] G_adj[i];
        delete [] G_adj;
    }

    if (used) {
        for (i = 0; i <= K; i++) delete [] used[i];
        delete [] used;
    }


    delete [] v_lab;

    delete [] e_lab;

    if (out_v_size) {
        for (i = 0; i <= K; i++) delete [] out_v_size[i];
        delete [] out_v_size;
    }

    if (out_e_size) {
        for (i = 0; i <= K; i++) delete [] out_e_size[i];
        delete [] out_e_size;
    }

    delete [] F;

    delete [] P;

    delete [] lack_size;

    if (lack) {
        for (i = 0; i < v_size; i++) delete [] lack[i];
        delete [] lack;
    }

    delete [] lev;

    delete [] loc;
}

void EBBkC_Graph_t::read_edges_from_file(const char *file_name) { /*该函数接受一个指向字符的指针作为参数，这个指针指向的是要读取的文件的名字。*/
    FILE *fp; /*声明了一个指向FILE的指针fp。*/
    if ((fp = fopen(file_name, "r")) == nullptr) { /*这行代码尝试以只读模式（"r"）打开名为file_name的文件，并将返回的文件指针赋给fp。*/
        printf("Cannot open file %s.\n", file_name); /*如果文件打开失败，这行代码会打印一条错误消息，消息内容包括无法打开的文件名。*/
        exit(0); /*终止程序的执行。*/
    }

    Edge_t e; /*声明了一个类型为Edge_t的变量e。*/
    int u, v, i;
    int *old2new = new int [N_NODES]; /*这个数组用于映射旧的节点编号到新的节点编号。*/
    for (i = 0; i < N_NODES; i++) old2new[i] = -1; /*它将old2new数组的所有元素初始化为-1。这意味着在一开始，所有的节点都还没有被映射到新的编号。*/

    e_size = 0; /*用于记录已经读取并处理的边的数量。*/
    edges = new Edge_t [N_EDGES];  /*用于存储读取并处理后的边。*/

    while (fscanf(fp, "%d %d%*[^\n]%*c", &u, &v) == 2) { /*它尝试从文件中读取两个整数，并将它们分别赋给变量u和v。如果成功读取两个整数，
    则继续执行循环体内的代码；否则退出循环。其中"%d %d"表示读取两个整数，"%[^\n]"表示跳过一行中的剩余部分（不包括换行符），"%*c"表示读取并丢弃换行符。
    这样可以处理文件中每行包含多个数字，但只有前两个数字被认为是边的情况。*/

        if (u > N_NODES || v > N_NODES) {/*这行代码检查读取到的节点编号是否超出了允许的范围。如果任何一个节点编号大于了预定义的常量N_NODES，则执行下面的代码块。*/
            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v)); /*打印一条消息，提示用户需要增大N_NODES的值，至少要大到文件中的最大节点编号。*/
            exit(0); /*如果节点编号超出了范围，这行代码会终止程序的执行。*/
        }
        if (old2new[u] == -1) { /*检查节点u是否已经被映射到一个新的编号。如果old2new[u]的值为-1，表示节点u还没有被映射。*/
            old2new[u] = (int) new2old.size(); /*如果节点u还没有被映射，这行代码将u映射到一个新的编号，这个新的编号就是new2old数组的大小。*/
            new2old.push_back(u); /*将节点u添加到new2old数组的末尾。我们可以推测new2old数组用于存储已经映射到新编号的节点。*/
        }
        if (old2new[v] == -1) { /*接下来的几行代码与上面的代码类似，但是它们是处理节点v而不是节点u。*/
            old2new[v] = (int) new2old.size();
            new2old.push_back(v);
        }

        e = Edge_t(old2new[u], old2new[v], false); /*创建了一个新的边对象e，它的两个节点分别是节点u和节点v映射后的新编号，它的第三个参数是false。
        我们可以推测这个第三个参数可能表示这条边是否是一个有效的边或者是否已经被处理过。*/
        edges[e_size++] = e; /*将新创建的边对象e添加到edges数组的末尾，并将e_size的值加1。这样，edges数组就存储了所有已经读取并处理的边，
        而e_size则记录了已经处理的边的数量。*/
    }

    v_size = (int) new2old.size(); /*计算并存储了已经映射到新编号的节点的数量。*/

    fclose(fp); /*关闭了之前打开的文件。*/

    delete [] old2new; /*释放了之前分配的用于存储节点映射关系的数组的内存空间。*/
}

void EBBkC_Graph_t::read_ordered_edges_from_file(const char *file_name) {/*这个函数的目的是从文件中读取有序的边。*/
    FILE *fp;/*定义一个文件指针fp。*/
    if ((fp = fopen(file_name, "r")) == nullptr) { /*尝试以只读模式打开指定的文件。如果文件打开失败，fopen函数返回nullptr。*/
        printf("Cannot open file %s.\n", file_name); /*如果文件打开失败，打印错误消息。*/
        exit(0); /*退出程序。*/
    }

    Edge_t e, e_; /*定义两个Edge_t类型的变量，可能表示图的边。*/
    int u, v, w, i, j, k, idx, edge_rank, edge_sub_size; /*定义多个整数变量，用于后续的读取和计算。*/
    int *old2new = new int [N_NODES]; /*动态分配一个整数数组，用于映射旧的节点编号到新的节点编号。*/
    for (i = 0; i < N_NODES; i++) old2new[i] = -1; /*初始化old2new数组，所有元素设置为-1。*/
    vector<int> t_; /*定义两个向量：一个整数向量和一个整数向量的向量。它们可能用于存储子图或子边的信息。*/
    vector<vector<int>> T_;

    e_size = 0; /*初始化边的大小为0。*/
    edges = new Edge_t [N_EDGES]; /*动态分配一个Edge_t数组，用于存储所有的边。*/

    while (fscanf(fp, "%d %d %d %d %d", &u, &v, &k, &edge_rank, &edge_sub_size) == 5) { /*它从文件中读取数据，直到文件结束或数据格式不符合预期：*/
    /*使用fscanf函数从文件中读取5个整数。如果成功读取5个整数，则继续循环。否则，退出循环。*/
        if (k <= K) { /*如果k小于或等于K，则跳过当前行的其余部分并继续下一次循环。这里假设我们只对k大于K的边感兴趣。*/
            fscanf(fp, "%*[^\n]%*c");
            continue;
        }
        truss_num = truss_num > k ? truss_num : k; /*使用三元运算符来确定truss_num的值：如果当前的k值大于truss_num，则更新truss_num为k的值；否则，保持原值。*/

        if (u > N_NODES || v > N_NODES) { /*检查读取到的节点编号u和v是否超出了预期的节点数量：*/
            printf("Enlarge N_NODES to at least %u.\n", (u > v ? u : v)); /*如果u或v大于N_NODES，则打印错误消息并退出程序。*/
            exit(0);
        }
        /*下面的代码段处理节点映射.对于节点u和v，如果它们在old2new数组中的值为-1（表示还没有映射），则为它们分配一个新的映射，并将它们添加到new2old向量中。*/
        if (old2new[u] == -1) {
            old2new[u] = (int) new2old.size();
            new2old.push_back(u);
        }
        if (old2new[v] == -1) {
            old2new[v] = (int) new2old.size();
            new2old.push_back(v);
        }
        /*接下来，代码创建一个新的边，并将其添加到数据结构中。使用映射后的节点编号创建一个新的边e，并将其添加到edges数组中。同时，更新边的ID映射和边的排名。*/
        e = Edge_t(old2new[u], old2new[v], false); /*这行代码创建了一个新的Edge_t对象e。这个对象接受三个参数：old2new[u]、old2new[v]和false。
        根据之前的代码，old2new是一个映射数组，它将旧的节点编号映射到新的节点编号。因此，old2new[u]和old2new[v]分别是节点u和v的新编号。*/
        rank.push_back(edge_rank); /*将edge_rank添加到rank向量的末尾。*/
        edge2id.insert(e, e_size); /*在edge2id映射中插入一个新的键值对。键是刚刚创建的边e，而值是e_size，它可能表示边在edges数组中的索引或位置。*/
        edges[e_size++] = e; /*它将新创建的边e添加到edges数组的当前位置（由e_size指定）。然后，它增加e_size的值，为下一条边在数组中的位置做准备。*/
        /*最后，代码处理与当前边相关的子边：*/
        for (j = 0; j < edge_sub_size; j++) { /*每次迭代处理一个子边。*/
            fscanf(fp, "%d", &w); /*这行代码从文件中读取一个整数，并将其存储在变量w中。*/
            if (old2new[w] == -1) { /*检查节点w是否已经被映射（即它在old2new数组中的值是否为-1）。如果没有映射，那么代码将为它分配一个新的映射，
            并将它添加到new2old向量中。*/
                old2new[w] = (int) new2old.size();
                new2old.push_back(w);
            }
            t_.push_back(old2new[w]); /*将节点w的新编号添加到t_向量的末尾。这个向量可能用于临时存储与当前边相关的子边的节点编号。*/
        }
        T_.push_back(t_); /*将整个t_向量（包含所有子边的节点编号）添加到T_向量的末尾。这意味着我们已经完成了当前边所有子边的处理。*/
        t_.clear(); /*清除t_向量的内容，为处理下一条边的子边做准备。*/
    }

    v_size = (int) new2old.size(); /*获取new2old向量的大小，并将其存储在v_size变量中。这代表了图中节点的数量。*/

    T = new int* [e_size]; /*用于存储与边相关的子边的节点编号。*/
    for (i = 0; i < e_size; i++) T[i] = new int [truss_num + 1]; /*为T数组中的每个元素动态分配一个整数数组，其大小为truss_num + 1。这个内部数组将用于存储特定边的子边的节点编号。*/
    T_size = new int [e_size]; /*这个数组将用于存储与每条边相关的子边的数量。*/

    C = new int* [e_size]; /*用于存储与每条边相关的候选边的索引。*/
    for (i = 0; i < e_size; i++) C[i] = new int [T_[i].size() * (T_[i].size() - 1) / 2];  /*这个循环为C数组中的每个元素动态分配一个整数数组，
    其大小是基于与边相关的子边的数量的组合数。这个内部数组将用于存储特定边的候选边的索引。*/
    C_size = new int [e_size]; /*用于存储与每条边相关的候选边的数量。*/

    for (i = 0; i < e_size; i++) { /*遍历所有的边，*/
        T_size[i] = C_size[i] = 0; /*初始化与当前边相关的子边和候选边的数量为0。*/

        for (j = 0; j < T_[i].size(); j++) T[i][T_size[i]++] = T_[i][j]; /*将与当前边相关的子边的节点编号从T_复制到T，并更新子边的数量。*/

        for (j = 0; j < T_size[i]; j++) { /*这个内部循环遍历与当前边相关的子边，*/
            for (k = j + 1; k < T_size[i]; k++) { /*循环遍历与当前子边配对的另一条子边*/

                e = Edge_t(T[i][j], T[i][k], false); /*创建一个新的边对象，其节点为当前子边和配对子边的节点。*/

                if ((idx = edge2id.exist(e)) != -1  && rank[idx] > rank[i]) { /*检查新的边对象是否在edge2id映射中存在，并且其排名是否高于当前边的排名。*/
                    C[i][C_size[i]++] = idx; /*如果条件满足，则将新的边的索引添加到当前边的候选边列表中，并更新候选边的数量。*/
                }
            }
        }
    }

    printf("|V| = %d, |E| = %d\n", v_size, e_size); /*打印图中节点和边的数量。这里使用了|V|和|E|来表示节点和边的数量。*/
    printf("Truss number = %d\n", truss_num - 2); /*打印Truss数，但是减去了2。*/

    fclose(fp); /*关闭文件指针，结束文件的读取操作。*/
    delete [] old2new; /*释放动态分配的old2new数组的内存空间。*/
}

void EBBkC_Graph_t::truss_decompose(const char* w_file_name) { /*该函数用于图的k-truss分解。*/
 
    int i, j, k, s, t, w, end, edge_seq = 1, sw, wt; /*定义了一系列整型变量，用于后续的循环和计算。*/

    auto *_d = new int [v_size](); /*这个数组用于存储每个节点的度数。*/
    auto *_cd = new int [v_size + 1]; /*这个数组用于存储每个节点的累积度数。*/
    auto *_adj = new int [2 * e_size]; /*这个数组用于存储邻接信息。*/

    for (i = 0; i < e_size; i++) { /*遍历每条边，更新每个节点的度数。*/
        _d[edges[i].s]++;
        _d[edges[i].t]++;
    }
    _cd[0] = 0;
    for (i = 1; i < v_size + 1; i++) { /*计算并存储每个节点的累积度数。同时将度数数组清零。*/
        _cd[i] = _cd[i - 1] + _d[i - 1];
        _d[i - 1] = 0;
    }
    for (i = 0; i < e_size; i++) { /*根据累积度数和度数信息，构建邻接数组。*/
        _adj[_cd[edges[i].s] + _d[edges[i].s]++] = edges[i].t;
        _adj[_cd[edges[i].t] + _d[edges[i].t]++] = edges[i].s;
    }

    KeyVal_t kv;
    Heap_t heap;
    Edge_t  e;

    auto *sup = new int [e_size](); /*动态分配了两个数组，一个用于存储支持度信息，另一个用于存储布尔值，可能用于检查节点是否已经被访问过。*/
    auto *h_table = new bool [v_size]();
    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_id; /*定义了一个无序映射，用于存储边的ID信息。*/

    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_truss; /*定义了三个无序映射，分别用于存储边的truss值、rank值和子图信息。*/
    unordered_map<Edge_t, int, Edge_t::Hash_Edge_t> edge_rank;
    unordered_map<Edge_t, vector<int>, Edge_t::Hash_Edge_t> edge_sub;

    for (i = 0; i < e_size; i++) { /*遍历每条边并执行一系列操作。具体操作包括：设置边的ID、根据节点度数选择s和t、通过邻接数组设置h_table、计算支持度、
    并重置h_table。这个循环的目的是计算每条边的支持度，这是k-truss分解的关键步骤之一。*/

        edge_id[edges[i]] = i; /*将每条边的ID存储到edge_id映射中。*/

        s = _d[edges[i].s] < _d[edges[i].t] ? edges[i].s : edges[i].t; /*选择度数较小的节点作为s，度数较大的节点作为t。如果度数相同，则任意选择。*/
        t = _d[edges[i].s] < _d[edges[i].t] ? edges[i].t : edges[i].s;

        for (j = _cd[s]; j < _cd[s + 1]; j++) { /*通过邻接数组，遍历节点s的所有邻居，并将它们标记为已访问。*/
            w = _adj[j];
            h_table[w] = true;
        }

        for (j = _cd[t]; j < _cd[t + 1]; j++) { /*遍历节点t的所有邻居，如果某个邻居已经被标记为已访问（即在s的邻居中），则增加边i的支持度。*/
            w = _adj[j];
            if (h_table[w])
                sup[i]++;
        }

        for (j = _cd[s]; j < _cd[s + 1]; j++) {/*重置h_table，将节点s的所有邻居标记为未访问。*/
            w = _adj[j];
            h_table[w] = false;
        }
    }

    heap.make_heap(sup, e_size); /*使用给定的支持度数组sup和边数量e_size来构建一个最小堆。*/

    for (k = 3; !heap.empty(); k++) { /*从k=3开始，直到堆为空为止。每次迭代处理一个特定的k-truss。*/

        while (!heap.empty() && heap.min_element().val < k - 2) { /*循环直到堆为空或者堆的最小元素的值不小于k-2。这是为了确保处理的是k-truss。*/
            kv = heap.pop(); /*从堆中弹出最小元素。*/
            e = edges[kv.key]; /*获取与弹出的键值对应的边。*/
            edge_truss[e] = k; /*设置当前边的truss值为k。*/
            edge_rank[e] = edge_seq++; /*为当前边分配一个唯一的序列号。*/

            s = _d[e.s] < _d[e.t] ? e.s : e.t; /*选择度数较小的节点作为s，度数较大的节点作为t。这部分与之前的代码段相似。*/
            t = _d[e.s] < _d[e.t] ? e.t : e.s;

            end = _cd[s] + _d[s]; /*计算节点s的邻接列表的结束位置。*/
            for (j = _cd[s]; j < end; j++) { /*遍历节点s的所有邻居，并将它们标记为已访问。*/
                w = _adj[j];
                h_table[w] = true;
            }

            end = _cd[t] + _d[t]; /*计算节点t的邻接列表的结束位置。*/
            for (j = _cd[t]; j < end; j++) { /*遍历节点t的所有邻居。*/
                w = _adj[j]; /*从邻接列表中获取当前邻居节点的ID，并将其存储在变量w中。*/

                if (w == s) { /*这个if语句检查当前邻居节点是否是节点s。如果是，执行以下操作：*/
                    _adj[j--] = _adj[--end]; /*将邻接列表中当前位置的元素替换为列表末尾的元素，并同时更新j和end的值。这是一种常见的在遍历过程中删除元素的技巧。*/
                    _d[t]--; /*减少节点t的度数，因为已经删除了一个邻居节点。*/
                }

                if (h_table[w]) { /*检查当前邻居节点w是否已被访问（在之前的代码中，如果节点被访问，其对应的h_table值会被设置为true）。如果已被访问，执行以下操作*/
                    edge_sub[e].push_back(w); /*将当前邻居节点添加到边e的子图列表中。*/
                    heap.update(edge_id[Edge_t(s, w, false)]); /*更新堆中与边(s, w)相关的元素的支持度。*/
                    heap.update(edge_id[Edge_t(w, t, false)]); /*更新堆中与边(w, t)相关的元素的支持度。*/
                }
            }

            end = _cd[s] + _d[s]; /*这部分代码与之前的类似，但是针对的是节点s的邻接列表。它重置了所有已访问邻居节点的状态，并在必要时更新节点s的邻接列表和度数。*/
            for (j = _cd[s]; j < end; j++) {
                w = _adj[j];
                h_table[w] = false; /*将已访问的邻居节点标记为未访问状态。*/

                if (w == t) { /*与之前的if语句类似，这个语句检查当前邻居节点是否是节点t，并在必要时更新邻接列表和度数。*/
                    _adj[j--] = _adj[--end];
                    _d[s]--;
                }
            }
        }
    }

    heap.release_heap(); /*释放堆的内存。*/

    FILE *fp = fopen(w_file_name, "w"); /*打开一个文件，准备写入结果。*/

    for (i = 0; i < e_size; i++) { /*遍历每条边，并将结果写入文件。结果包括边的两个节点、truss值、rank值和子图信息。*/
        e = edges[i];
        fprintf(fp, "%d %d %d %d %ld ", new2old[e.s], new2old[e.t], edge_truss[e], edge_rank[e], edge_sub[e].size());
        for (j = 0; j < edge_sub[e].size(); j++) fprintf(fp, "%d ", new2old[edge_sub[e][j]]);
        fprintf(fp, "\n");
    }

    fclose(fp); /*关闭文件。*/

    delete [] _d; /*删除动态分配的内存空间。*/
    delete [] _cd;
    delete [] _adj;
    delete [] h_table;
}


void EBBkC_Graph_t::build_from_G() {
    int i;

    sub_v = new int* [K + 1];

    sub_e = new int* [K + 1];

    sub_e_size = new int [K + 1];

    sub_v_size = new int [K + 1];

    for (i = 0; i < K; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[K] = new int [v_size];

    for (i = 0; i < K; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[K] = new int [e_size];

    sub_v_size[K] = 0;
    for (i = 0; i < v_size; i++) sub_v[K][sub_v_size[K]++] = i;

    sub_e_size[K] = 0;
    for (i = 0; i < e_size; i++) sub_e[K][sub_e_size[K]++] = i;

    lab = new int [v_size];
    for (i = 0; i < v_size; i++) lab[i] = K;

    DAG_deg = new int* [K + 1];
    for (i = 0; i <= K; i++) DAG_deg[i] = new int [v_size];

    G_deg = new int* [K + 1];
    for (i = 0; i <= K; i++) G_deg[i] = new int [v_size];

    col = new int [v_size];

    DAG_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) DAG_adj[i] = new int [truss_num + 1];

    G_adj = new int* [v_size];
    for (i = 0; i < v_size; i++) G_adj[i] = new int  [truss_num + 1];

    used = new bool* [K + 1];
    for (i = 0; i <= K; i++) used[i] = new bool [v_size + 1]();

    v_lab = new int [v_size];
    for (i = 0; i < v_size; i++) v_lab[i] = K;

    e_lab = new int [e_size];
    for (i = 0; i < e_size; i++) e_lab[i] = K;

    out_v_size = new int* [K + 1];
    for (i = 0; i <= K; i++) out_v_size[i] = new int [e_size];

    out_e_size = new int* [K + 1];
    for (i = 0; i <= K; i++) out_e_size[i] = new int [e_size];

    F = new int [truss_num + 1];

    P = new int [truss_num + 1];

    lack_size = new int [v_size];

    lack = new int* [v_size];
    for (i = 0; i < v_size; i++) lack[i] = new int [L + 1];

    lev = new int [v_size]();

    loc = new int [v_size];
}


void EBBkC_Graph_t::EBBkC(int l, unsigned long long *cliques) {
    int i, j, k, u, e, e_, _e, end;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else if (K == 4) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else {

            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                sort(T[e], T[e] + T_size[e]);
                sort(C[e], C[e] + C_size[e]);
            }
            
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;

                sub_v_size[l - 2] = 0;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
                }

                sub_e_size[l - 2] = 0;

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    sub_e[l - 2][sub_e_size[l - 2]++] = e_;
                }

                EBBkC(l - 2, cliques);

            }
        }

        return;
    }

    for (i = 0; i < sub_e_size[l]; i++) {
        e = sub_e[l][i];

        if (l == 3) {
            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
            for (j = 0; j < sub_v_size[l - 2]; j++) {
                (*cliques)++;
            }
        }

        else if (l == 4) {
             sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
             for (j = 0; j < sub_e_size[l - 2]; j++) {
                (*cliques)++;
             }
        }

        else {
            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
            sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
            EBBkC(l - 2, cliques);
        }
    }
}


void EBBkC_Graph_t::EBBkC_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                 e = sub_e[l][i];

                 for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                 }
            }
        }
        else if (K == 4) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    col[u] = 0;
                    DAG_deg[0][u] = 0;
                    G_deg[l - 2][u] = 0;
                }

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    s = edges[e_].s;
                    t = edges[e_].t;
                    G_adj[s][G_deg[l - 2][s]++] = t;
                    G_adj[t][G_deg[l - 2][t]++] = s;
                }

                auto *list = new KeyVal_t [truss_num + 1];
                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    list[j].key = u;
                    list[j].val = G_deg[l - 2][u];
                }
                sort(list, list + T_size[e]);

                for (j = 0; j < T_size[e]; j++) {
                    u = list[j].key;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = true;
                    }
                    for (c = 1; used[K][c]; c++) ;
                    col[u] = c;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = false;
                    }
                }
                delete [] list;

                sub_v_size[l - 2] = 0;
                dist = 0;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
                    if (!used[K][col[u]]) {
                        used[K][col[u]] = true;
                        dist++;
                    }
                }

                if (dist >= l - 2) {
                    sort(sub_v[l - 2], sub_v[l - 2] + sub_v_size[l - 2]);

                    sub_e_size[l - 2] = 0;
                    for (j = 0; j < C_size[e]; j++) {
                        e_ = C[e][j];
                        sub_e[l - 2][sub_e_size[l - 2]++] = e_;
                        s = edges[e_].s;
                        t = edges[e_].t;
                        edges[e_].s = (col[s] > col[t]) ? s : t;
                        edges[e_].t = (col[s] > col[t]) ? t : s;
                        s = edges[e_].s;
                        t = edges[e_].t;
                        DAG_adj[s][DAG_deg[0][s]++] = t;
                    }

                    for (j = 0; j < T_size[i]; j++) {
                        u = T[e][j];
                        // sorted array for SIMD usage, DAG_adj can be considered as const.
                        sort(DAG_adj[u], DAG_adj[u] + DAG_deg[0][u]);
                    }

                    EBBkC_plus(l - 2, cliques);
                }

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    used[K][col[u]] = false;
                }
            }
        }

        return;
    }

    if (l == 1) {
        for (i = 0; i < sub_v_size[l]; i++) {
            (*cliques)++;
        }

        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = intersect_simd4x(sub_v[l], sub_v_size[l], DAG_adj[u], DAG_deg[0][u], sub_v[l - 1]);

        if (l == 2) {
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                (*cliques)++;
            }
        }

        else {
            if (sub_v_size[l - 1] >= l - 1) {
                EBBkC_plus(l - 1, cliques);
            }
        }
    }
}


void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist;

    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) {
        if (K == 3) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];
                for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else if (K == 4) {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];
                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else {
            for (i = 0; i < sub_e_size[l]; i++) {
                e = sub_e[l][i];

                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    col[u] = 0;
                    DAG_deg[l - 2][u] = 0;
                    G_deg[l - 2][u] = 0;
                }

                for (j = 0; j < C_size[e]; j++) {
                    e_ = C[e][j];
                    s = edges[e_].s;
                    t = edges[e_].t;
                    G_adj[s][G_deg[l - 2][s]++] = t;
                    G_adj[t][G_deg[l - 2][t]++] = s;
                }

                auto *list = new KeyVal_t [truss_num + 1];
                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    list[j].key = u;
                    list[j].val = G_deg[l - 2][u];
                }
                sort(list, list + T_size[e]);

                for (j = 0; j < T_size[e]; j++) {
                    u = list[j].key;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = true;
                    }
                    for (c = 1; used[K][c]; c++) ;
                    col[u] = c;
                    for (k = 0; k < G_deg[l - 2][u]; k++) {
                        v = G_adj[u][k];
                        used[K][col[v]] = false;
                    }
                }
                delete [] list;

                sub_v_size[l - 2] = 0;
                dist = 0;

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u;
                    if (!used[K][col[u]]) {
                        used[K][col[u]] = true;
                        dist++;
                    }
                }

                if (dist >= l - 2) {
                    sub_e_size[l - 2] = 0;
                    for (j = 0; j < C_size[e]; j++) {
                        e_ = C[e][j];
                        sub_e[l - 2][sub_e_size[l - 2]++] = e_;
                        s = edges[e_].s;
                        t = edges[e_].t;
                        edges[e_].s = (col[s] > col[t]) ? s : t;
                        edges[e_].t = (col[s] > col[t]) ? t : s;
                        s = edges[e_].s;
                        t = edges[e_].t;

                        DAG_adj[s][DAG_deg[l - 2][s]++] = t;
                    }

                    EBBkC_plus_plus(l - 2, cliques);
                }

                for (j = 0; j < T_size[e]; j++) {
                    u = T[e][j];
                    used[K][col[u]] = false;
                }
            }
        }

        return;
    }

    if (l == 2) {
        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];

            for (j = 0; j < DAG_deg[l][u]; j++) {
                (*cliques)++;
            }
        }

        return;
    }

    if (l == 3) {

        for (i = 0; i < sub_v_size[l]; i++) {
            u = sub_v[l][i];

            if (col[u] < l) continue;

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                lab[v] = l - 1;
            }

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];

                if (col[v] < l - 1) continue;

                for (k = 0; k < DAG_deg[l][v]; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) (*cliques)++;
                }
            }

            for (j = 0; j < DAG_deg[l][u]; j++) {
                v = DAG_adj[u][j];
                lab[v] = l;
            }
        }

        return;
    }

    if (can_terminate(l, cliques)) {
        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (col[u] < l) continue;

        sub_v_size[l - 1] = 0;
        dist = 0;

        for (j = 0; j < DAG_deg[l][u]; j++) {
            v = DAG_adj[u][j];
            lab[v] = l - 1;
            sub_v[l - 1][sub_v_size[l - 1]++] = v;
            DAG_deg[l - 1][v] = 0;
            G_deg[l - 1][v] = 0;

            if (!used[l][col[v]]) {
                used[l][col[v]] = true;
                dist++;
            }
        }

        if (dist >= l - 1) {

            sub_e_size[l - 1] = 0;
            for (j = 0; j < sub_v_size[l - 1]; j++) {
                v = sub_v[l - 1][j];

                end = DAG_deg[l][v];
                for (k = 0; k < end; k++) {
                    w = DAG_adj[v][k];
                    if (lab[w] == l - 1) {
                        DAG_deg[l - 1][v]++;
                        sub_e_size[l - 1]++;

                        // just for early-termination
                        G_deg[l - 1][v]++;
                        G_deg[l - 1][w]++;

                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end];
                        DAG_adj[v][end] = w;
                    }
                }
            }

            EBBkC_plus_plus(l - 1, cliques);
        }

        for (j = 0; j < sub_v_size[l - 1]; j++) {
            v = sub_v[l - 1][j];
            lab[v] = l;
            used[l][col[v]] = false;
        }
    }
}


void EBBkC_Comb_list(int *list, int list_size, int start, int picked, int k, unsigned long long *cliques) {
    if (picked == k) {
        (*cliques)++;
        return;
    }

    for (int i = start; i < list_size; i++) {
        EBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques);
    }
}

void EBBkC_Graph_t::list_in_plex(int start, int p, int q, unsigned long long *cliques) {
    if (F_size < q) return;

    if (p == 0) {
        if (q > F_size - q)
            EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques);
        else
            EBBkC_Comb_list(F, F_size, 0, 0, q, cliques);
        return;
    }

    int i, j, u, v, vis = 0;

    for (i = start; i < P_size && P_act >= p; i++) {
        u = P[i];

        if (lev[u]) continue;

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == 0) {
                lev[v] = p;
                P_act--;
            }
        }

        list_in_plex(i + 1, p - 1, q, cliques);

        for (j = 0; j < lack_size[u]; j++) {
            v = lack[u][j];
            if (loc[v] >= i && lev[v] == p) {
                lev[v] = 0;
                P_act++;
            }
        }

        P_act--;
        vis++;
    }
    P_act += vis;
}

bool EBBkC_Graph_t::can_terminate(int l, unsigned long long *cliques) {
    int i, j, k, u, v, end, p_;

    if (sub_e_size[l] < sub_v_size[l] * (sub_v_size[l] - L) / 2) return false;

    if (sub_e_size[l] == sub_v_size[l] * (sub_v_size[l] - 1) / 2) {
        if (l > sub_v_size[l] - l)
            EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, sub_v_size[l] - l, cliques);
        else
            EBBkC_Comb_list(sub_v[l], sub_v_size[l], 0, 0, l, cliques);

        return true;
    }

    if (L == 1) return false;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (sub_v_size[l] - G_deg[l][u] > L) {
            return false;
        }
    }

    F_size = 0;
    P_size = 0;

    for (i = 0; i < sub_v_size[l]; i++) {
        u = sub_v[l][i];

        if (G_deg[l][u] == sub_v_size[l] - 1) {
            loc[u] = -1;
            F[F_size++] = u;
            continue;
        }

        loc[u] = P_size;
        P[P_size++] = u;
    }

    int* e = new int [P_size * P_size]();

    for (i = 0; i < P_size; i++) {
        u = P[i];
        lack_size[u] = 0;

        end = DAG_deg[l][u];
        for (j = 0; j < end; j++) {
            v = DAG_adj[u][j];

            if (loc[v] != -1) {
                e[loc[u] * P_size + loc[v]] = 1;
                e[loc[v] * P_size + loc[u]] = 1;
            }
        }
    }

    for (i = 0; i < P_size * P_size; i++) {
        if (!e[i]) {
            j = i / P_size, k = i % P_size;
            u = P[j], v = P[k];
            lack[u][lack_size[u]++] = v;
        }
    }

    delete [] e;

    for(i = 0; i <= l; i++) {
        P_act = P_size;
        list_in_plex(0, i, l - i, cliques);
    }

    return true;
}


double EBBkC_t::truss_order(const char *r_file_name, const char *w_file_name) { /*接受两个参数：输入文件名r_file_name和输出文件名w_file_name。返回类型为double。*/
    double runtime; /*用于存储运行时间。*/
    struct rusage start, end; /*这两个结构体通常用于存储程序运行的时间消耗。*/
    EBBkC_Graph_t G; /*用于存储图的有关信息。*/

    GetCurTime(&start); /*获取当前时间并存储在变量start中。*/
    G.read_edges_from_file(r_file_name); /*从给定的文件名读取边的信息并存储在G中。*/
    G.truss_decompose(w_file_name); /*对G进行Truss分解，并将结果（可能是一个子图或者边的某种分解）输出到给定的文件名。*/
    GetCurTime(&end); /*再次调用GetCurTime函数，获取当前时间并存储在变量end中。这一步是为了获取执行truss_decompose函数后的时间。*/
    runtime = GetTime(&start, &end);/*获取start和end之间的时间差，并将结果存储在runtime变量中。这个时间差就是执行truss_decompose函数的运行时间。*/

    return runtime; /*返回运行时间runtime。*/
}

double EBBkC_t::list_k_clique(const char *file_name) { /*接受一个文件名字符串作为参数，返回一个 double 类型的值，代表运行时间。*/
    double runtime;  /*用于存储函数的运行时间。*/
    struct rusage start, end; /*定义了两个 rusage 结构体变量 start 和 end，用于存储函数开始和结束时的资源使用情况。*/
    EBBkC_Graph_t G; /*用于表示图的自定义数据结构。*/

    printf("Reading edges from %s ...\n", file_name); /*打印消息到控制台，表示正在从指定的文件名读取边。*/
    G.read_ordered_edges_from_file(file_name); /*从文件中读取有序的边。传入的文件名是从函数的参数获得的。*/

    printf("Building necessary data structure ...\n"); /*打印消息到控制台，表示正在构建必要的数据结构。*/
    G.build_from_G(); /*根据已读取的边构建必要的数据结构。*/

    printf("Iterate over all cliques\n"); /*打印消息到控制台，表示将遍历所有的k-clique。*/

    GetCurTime(&start); /*获取当前的时间使用情况，并将其存储在 start 变量中。*/
    G.EBBkC_plus_plus(K, &N); /*执行k-clique的查找算法。*/
    GetCurTime(&end); /*再次获取当前的时间使用情况，并将其存储在 end 变量中。*/
    runtime = GetTime(&start, &end); /*执行k-clique查找算法所花费的时间。*/

    return runtime; /*返回函数的运行时间。*/
}