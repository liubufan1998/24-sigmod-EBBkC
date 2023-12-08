#include "edge_oriented.h"
#include <set>
#include <algorithm>
#include <unordered_map>
#include "set_operation.h"

extern const int K, L; /*这两行声明了外部常量K、L和N，这些常量可能在类的其他部分或全局范围内定义。*/
extern unsigned long long N;

EBBkC_Graph_t::EBBkC_Graph_t() = default; /*这是类的默认构造函数的定义，使用默认关键字表示使用编译器生成的默认构造函数。*/

EBBkC_Graph_t::~EBBkC_Graph_t() { /*类的析构函数的开始，用于释放对象占用的资源。*/
    int i;

    delete [] edges; /*释放edges数组占用的内存。*/

    if (T) { /*检查成员变量T是否为非空，如果是，则遍历并删除每个子数组，然后删除数组T本身。*/
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


void EBBkC_Graph_t::build_from_G() { /*它的主要目的是为一个图数据结构分配必要的内存，并初始化一些变量。*/
    int i; /*声明一个整型变量i，用于后续的循环。*/

    sub_v = new int* [K + 1]; /*为sub_v和sub_e分配内存，它们是指向指针的指针，用于存储子图和子边的信息。*/

    sub_e = new int* [K + 1];

    sub_e_size = new int [K + 1]; /*它们分别存储子边和子顶点的数量。*/

    sub_v_size = new int [K + 1];
    /*使用循环为sub_v的前K个子数组分配内存。每个子数组的大小基于truss_num或v_size*/
    for (i = 0; i < K; i++) sub_v[i] = new int [truss_num + 1];
    sub_v[K] = new int [v_size];
    /*使用循环为sub_e的前K个子数组分配内存。每个子数组的大小基于truss_num或e_size。*/
    for (i = 0; i < K; i++) sub_e[i] = new int [truss_num * (truss_num - 1) / 2];
    sub_e[K] = new int [e_size];
    /*初始化sub_v_size[K]和sub_e_size[K]，并将所有的顶点和边添加到相应的子图和子边中。*/
    sub_v_size[K] = 0;
    for (i = 0; i < v_size; i++) sub_v[K][sub_v_size[K]++] = i;

    sub_e_size[K] = 0;
    for (i = 0; i < e_size; i++) sub_e[K][sub_e_size[K]++] = i;

    lab = new int [v_size]; /*为顶点标签数组分配内存。*/
    for (i = 0; i < v_size; i++) lab[i] = K; /*使用循环初始化所有顶点的标签为K。*/

    DAG_deg = new int* [K + 1]; /*为DAG（有向无环图）的度数数组分配内存。*/
    for (i = 0; i <= K; i++) DAG_deg[i] = new int [v_size]; /* 使用循环为DAG度数数组的每一个元素分配内存。*/

    G_deg = new int* [K + 1]; /*为G的度数数组分配内存。*/
    for (i = 0; i <= K; i++) G_deg[i] = new int [v_size]; /*使用循环为G度数数组的每一个元素分配内存。*/

    col = new int [v_size]; /*为颜色数组分配内存。*/
 
    DAG_adj = new int* [v_size]; /*为DAG的邻接数组分配内存。*/
    for (i = 0; i < v_size; i++) DAG_adj[i] = new int [truss_num + 1]; /*使用循环为DAG邻接数组的每一个元素分配内存。*/

    G_adj = new int* [v_size]; /*使用循环为G邻接数组的每一个元素分配内存。*/
    for (i = 0; i < v_size; i++) G_adj[i] = new int  [truss_num + 1]; /*为used数组分配内存，它是一个布尔型数组，用于跟踪哪些顶点或边已被使用。*/

    used = new bool* [K + 1]; /*使用循环初始化used数组的每一个元素。*/
    for (i = 0; i <= K; i++) used[i] = new bool [v_size + 1]();

    v_lab = new int [v_size]; /*为顶点标签数组和边标签数组分配内存，并使用循环进行初始化。*/
    for (i = 0; i < v_size; i++) v_lab[i] = K;

    e_lab = new int [e_size];
    for (i = 0; i < e_size; i++) e_lab[i] = K;

    out_v_size = new int* [K + 1];/*为输出顶点大小和输出边大小数组分配内存。*/
    for (i = 0; i <= K; i++) out_v_size[i] = new int [e_size];

    out_e_size = new int* [K + 1]; /*使用循环为输出顶点大小和输出边大小数组的每一个元素分配内存。*/
    for (i = 0; i <= K; i++) out_e_size[i] = new int [e_size];

    F = new int [truss_num + 1]; /*为F和P数组分配内存，它们可能是用于存储某种临时信息的数组。*/

    P = new int [truss_num + 1];

    lack_size = new int [v_size]; /*为lack_size数组分配内存，它可能用于跟踪每个顶点的某种缺乏的数量。*/

    lack = new int* [v_size]; /*为lack数组分配内存，它是一个整数指针数组，可能用于存储与每个顶点相关的某种缺乏的信息。*/
    for (i = 0; i < v_size; i++) lack[i] = new int [L + 1]; /*使用循环为lack数组的每一个元素分配内存。*/

    lev = new int [v_size](); /*为lev数组分配内存并使用括号初始化所有元素为0，它可能用于存储每个顶点的级别信息。*/

    loc = new int [v_size]; /*为loc数组分配内存，它可能用于存储每个顶点的位置信息。*/
}


void EBBkC_Graph_t::EBBkC(int l, unsigned long long *cliques) { /*接受两个参数：一个整数l和一个指向无符号长长整型的指针cliques。*/
    int i, j, k, u, e, e_, _e, end; /*声明了一些整型变量，用于后续的循环和临时存储。*/
    /*如果子图的顶点数量小于l或者子图的边数量小于l * (l - 1) / 2（即l个顶点可以形成的最大边数），则函数直接返回。*/
    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return;

    if (l == K) { /*再次判断，如果l等于K，则执行以下的代码块。*/
        if (K == 3) { /*如果K等于3，执行以下的代码块。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*外层循环遍历所有的子边，内层循环根据每条边的T_size来增加cliques的计数。*/
                e = sub_e[l][i];
                for (j = 0; j < T_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else if (K == 4) { /*如果K等于4，执行以下的代码块。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*与上面的代码块类似，但是这次是根据每条边的C_size来增加cliques的计数。*/
                e = sub_e[l][i];
                for (j = 0; j < C_size[e]; j++) {
                    (*cliques)++;
                }
            }
        }
        else { /*与上面的代码块类似，但是这次是根据每条边的C_size来增加cliques的计数。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*遍历所有子图l中的边。*/
                e = sub_e[l][i]; /*获取子图l中的第i条边，并将其存储在变量e中.*/
                sort(T[e], T[e] + T_size[e]); /* 对边e对应的T数组进行排序。这里排序的是从T[e]开始到T[e] + T_size[e]的部分。*/
                sort(C[e], C[e] + C_size[e]); /*对边e对应的C数组进行排序。这里排序的是从C[e]开始到C[e] + C_size[e]的部分。*/
            }
            
            for (i = 0; i < sub_e_size[l]; i++) { /*遍历所有子图l中的边。*/
                e = sub_e[l][i];
                /*如果边e的T数组大小小于l-2或C数组大小小于(l-2)*(l-3)/2，则跳过此次循环。这是为了确保有足够的顶点和边来形成一个更小的K-Clique。*/
                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;

                sub_v_size[l - 2] = 0; /*重置子图l-2的顶点数量。*/

                for (j = 0; j < T_size[e]; j++) { /*循环遍历边e的T数组，并将顶点添加到子图l-2中。*/
                    u = T[e][j];
                    sub_v[l - 2][sub_v_size[l - 2]++] = u; /*将顶点u添加到子图l-2的顶点数组中，并递增顶点数量*/
                }

                sub_e_size[l - 2] = 0; /*重置子图l-2的边数量。*/

                for (j = 0; j < C_size[e]; j++) { /*循环遍历边e的C数组，并将边添加到子图l-2中。*/
                    e_ = C[e][j]; 
                    sub_e[l - 2][sub_e_size[l - 2]++] = e_; /*将边e_添加到子图l-2的边数组中，并递增边的数量。*/
                }

                EBBkC(l - 2, cliques); /*对子图l-2递归调用EBBkC函数，以寻找更小的K-Clique。*/

            }
        }

        return;
    }

    for (i = 0; i < sub_e_size[l]; i++) {
        e = sub_e[l][i]; /*获取子图l中的第i条边，并将其存储在变量e中。*/

        if (l == 3) { /*判断当前处理的子图的大小是否为3。*/
            /*如果子图大小为3，使用函数intersect_simd4x计算两个顶点集合的交集，并返回交集的大小。*/
            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
            for (j = 0; j < sub_v_size[l - 2]; j++) { /*遍历上一步计算的交集结果。*/
                (*cliques)++; /*对于交集中的每个顶点，增加cliques的计数。这表示找到了一个新的3-clique。*/
            }
        }

        else if (l == 4) { /*判断当前处理的子图的大小是否为4。*/
            /*如果子图大小为4，使用函数intersect_simd4x计算两个边集合的交集，并返回交集的大小。*/
             sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
             for (j = 0; j < sub_e_size[l - 2]; j++) { /*遍历上一步计算的交集结果。*/
                (*cliques)++; /* 对于交集中的每条边，增加cliques的计数。这表示找到了一个新的4-clique。*/
             }
        }

        else { /* 如果子图的大小不是3或4，执行以下代码块。*/
            /*计算子图l的顶点与边e的T数组的交集，结果存储在sub_v[l - 2]中。*/
            sub_v_size[l - 2] = intersect_simd4x(sub_v[l], sub_v_size[l], T[e], T_size[e], sub_v[l - 2]);
            /*计算子图l的边与边e的C数组的交集，结果存储在sub_e[l - 2]中。*/
            sub_e_size[l - 2] = intersect_simd4x(sub_e[l], sub_e_size[l], C[e], C_size[e], sub_e[l - 2]);
            EBBkC(l - 2, cliques); /*对更小的子图（大小为l-2）递归调用函数EBBkC以继续寻找cliques。*/
        }
    }
}


void EBBkC_Graph_t::EBBkC_plus(int l, unsigned long long *cliques) { /*它接受两个参数：一个整数 l 和一个指向无符号长长整型的指针 cliques。*/
    int c, i, j, k, e, e_, u, v, w, s, t, end, dist; /*声明了一系列整数变量，用于在函数内部进行各种计算。*/

    if (sub_v_size[l] < l) return; /*如果子图l的顶点数小于l，则函数立即返回。*/

    if (l == K) { /*判断当前处理的子图大小是否等于某个常数K。*/
        if (K == 3) { /*如果K等于3，执行以下代码块。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*遍历子图l中的所有边。*/
                 e = sub_e[l][i]; /*获取子图l中的第i条边，并将其存储在变量e中。*/

                 for (j = 0; j < T_size[e]; j++) { /*遍历与边e相关的T数组的所有元素。*/
                    (*cliques)++; /*对于T数组中的每个元素，增加cliques的计数。*/
                 }
            }
        }
        else if (K == 4) { /*如果K等于4，执行以下代码块。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*遍历子图l中的所有边。*/
                e = sub_e[l][i]; /*获取子图l中的第i条边，并将其存储在变量e中。*/

                for (j = 0; j < C_size[e]; j++) { /*遍历与边e相关的C数组的所有元素。*/
                    (*cliques)++; /*对于C数组中的每个元素，增加cliques的计数。*/
                }
            }
        }
        else { /*如果K既不等于3也不等于4，执行以下代码块:*/
            for (i = 0; i < sub_e_size[l]; i++) {/*遍历子图l中的所有边。*/
                e = sub_e[l][i]; /*获取子图l中的第i条边，并将其存储在变量e中。*/
                /* 如果与边e相关的T数组的大小小于l-2，或者C数组的大小小于(l-2)*(l-3)/2，则跳过当前迭代。*/
                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue; 

                for (j = 0; j < T_size[e]; j++) { /*遍历与边e相关的T数组的所有元素。*/
                    u = T[e][j]; /*获取T数组中的第j个元素，并将其存储在变量u中。*/
                    col[u] = 0; /*将顶点u的col属性设置为0。*/
                    DAG_deg[0][u] = 0; /*将顶点u的DAG_deg属性（似乎是一个二维数组）的第0行设置为0。*/
                    G_deg[l - 2][u] = 0; /*将顶点u的G_deg属性的第l-2行设置为0。*/
                }

                for (j = 0; j < C_size[e]; j++) { /*这个循环遍历与边e相关的C数组的所有元素。*/
                    e_ = C[e][j]; /*获取C数组中的第j个元素，并将其存储在变量e_中。*/
                    s = edges[e_].s; /*获取边e_的起始顶点，并将其存储在变量s中. */
                    t = edges[e_].t; /*获取边e_的终止顶点，并将其存储在变量t中。*/
                    G_adj[s][G_deg[l - 2][s]++] = t; /*在邻接矩阵G_adj中，将顶点s和t连接起来，并增加s的度数。*/
                    G_adj[t][G_deg[l - 2][t]++] = s; /*在邻接矩阵G_adj中，将顶点t和s连接起来，并增加t的度数。*/
                }

                auto *list = new KeyVal_t [truss_num + 1]; /*动态分配一个KeyVal_t类型的数组，并将其指针存储在变量list中。*/
                for (j = 0; j < T_size[e]; j++) { /*循环遍历与边e相关的T数组的所有元素。*/
                    u = T[e][j]; /*获取T数组中的第j个元素，并将其存储在变量u中。*/
                    list[j].key = u; /*将u设置为list数组中第j个元素的键。*/
                    list[j].val = G_deg[l - 2][u]; /* 将u的度数设置为list数组中第j个元素的值。*/
                }
                sort(list, list + T_size[e]); /*对list数组进行排序，按照度数从小到大的顺序排列。*/

                for (j = 0; j < T_size[e]; j++) { /*这个循环遍历排序后的list数组的所有元素。*/
                    u = list[j].key; /*获取排序后list数组中的第j个元素的键，并将其存储在变量u中。*/
                    /*接下来的几个循环和判断语句是用于对顶点进行着色（即分配一个唯一的颜色）的过程，其中使用了贪心算法的思想。*/
                    for (k = 0; k < G_deg[l - 2][u]; k++) { /*这个内层循环遍历与顶点u相邻的所有顶点。*/
                        v = G_adj[u][k]; /*获取u的邻接顶点v。*/
                        used[K][col[v]] = true; /*将v的颜色标记为已使用。*/
                    }
                    for (c = 1; used[K][c]; c++) ; /*这个循环用于查找下一个可用的颜色。*/
                    col[u] = c; /*将顶点u的颜色设置为c。*/
                    for (k = 0; k < G_deg[l - 2][u]; k++) { /*内层循环再次遍历与顶点u相邻的所有顶点。*/
                        v = G_adj[u][k]; /*获取u的邻接顶点v。*/
                        used[K][col[v]] = false; /*将v的颜色标记为未使用。*/
                    }
                }
                delete [] list; /*释放之前动态分配的list数组的内存。*/

                sub_v_size[l - 2] = 0; /* 将sub_v_size数组中索引为l - 2的元素设置为0，表示子图中顶点的数量为0。*/
                dist = 0; /* 将变量dist设置为0，用于记录后续处理中着色顶点的数量。*/

                for (j = 0; j < T_size[e]; j++) { /*循环遍历与边e相关的T数组的所有元素。*/
                    u = T[e][j]; /*获取T数组中的第j个元素，并将其存储在变量u中。*/
                    sub_v[l - 2][sub_v_size[l - 2]++] = u; /*将顶点u添加到sub_v数组中，并增加sub_v_size[l - 2]的值。*/
                    if (!used[K][col[u]]) { /*判断顶点u的颜色是否已被使用。*/
                        used[K][col[u]] = true; /* 如果顶点u的颜色未被使用，则将其标记为已使用。*/
                        dist++; /*增加着色顶点的数量。*/
                    }
                }

                if (dist >= l - 2) { /*判断着色顶点的数量是否大于等于l - 2。*/
                    sort(sub_v[l - 2], sub_v[l - 2] + sub_v_size[l - 2]); /* 如果条件满足，则对sub_v数组进行排序。*/

                    sub_e_size[l - 2] = 0; /*将sub_e_size数组中索引为l - 2的元素设置为0，表示子图中边的数量为0。*/
                    for (j = 0; j < C_size[e]; j++) { /*循环遍历与边e相关的C数组的所有元素。*/
                        e_ = C[e][j];/*获取C数组中的第j个元素，并将其存储在变量e_中。*/
                        sub_e[l - 2][sub_e_size[l - 2]++] = e_; /*将边e_添加到sub_e数组中，并增加sub_e_size[l - 2]的值。*/
                        s = edges[e_].s; /* 获取边e_的起始顶点，并将其存储在变量s中。*/
                        t = edges[e_].t; /*获取边e_的终止顶点，并将其存储在变量t中。*/
                        edges[e_].s = (col[s] > col[t]) ? s : t; /*根据顶点的颜色对边的起始顶点和终止顶点进行排序。*/
                        edges[e_].t = (col[s] > col[t]) ? t : s; /* 根据顶点的颜色对边的起始顶点和终止顶点进行排序。*/
                        s = edges[e_].s; /*重新获取排序后的边的起始顶点，并将其存储在变量s中。*/
                        t = edges[e_].t; /*重新获取排序后的边的终止顶点，并将其存储在变量t中。*/
                        DAG_adj[s][DAG_deg[0][s]++] = t; /*将排序后的边添加到DAG的邻接矩阵中，并更新起始顶点的度数。*/
                    }

                    for (j = 0; j < T_size[i]; j++) { /*循环遍历与边i相关的T数组的所有元素。*/
                        u = T[e][j]; /*获取T数组中的第j个元素，并将其存储在变量u中。*/
                        // sorted array for SIMD usage, DAG_adj can be considered as const.
                        sort(DAG_adj[u], DAG_adj[u] + DAG_deg[0][u]); /*对DAG邻接矩阵中与顶点u相关的边进行排序。注释提到这是为了SIMD
                        （单指令多数据）的使用，说明排序后的数组可以提高某些处理效率。*/
                    }

                    EBBkC_plus(l - 2, cliques); /*调用EBBkC_plus函数，参数为l - 2和cliques。*/
                }

                for (j = 0; j < T_size[e]; j++) { /*循环遍历与边e相关的T数组的所有元素。*/
                    u = T[e][j]; /*获取T数组中的第j个元素，并将其存储在变量u中。*/
                    used[K][col[u]] = false; /*将used数组中与顶点u的颜色相关的元素设置为false，表示该颜色现在未被使用。*/
                }
            }
        }
        return;
    }

    if (l == 1) { /*判断变量l是否等于1。*/
        for (i = 0; i < sub_v_size[l]; i++) { /*如果l等于1，则执行这个循环，遍历子图顶点数组sub_v的所有元素。*/
            (*cliques)++; /*对cliques指针指向的值进行自增操作。*/
        }/*这段代码的目的是当l等于1时，对cliques指针指向的值进行自增操作，自增的次数为子图顶点数组sub_v的大小，然后函数返回。*/
        return;
    }

    for (i = 0; i < sub_v_size[l]; i++) { /*循环遍历子图顶点数组sub_v的所有元素。*/
        u = sub_v[l][i]; /*获取子图顶点数组sub_v的第i个元素，并将其存储在变量u中。*/

        if (col[u] < l) continue; /*判断顶点u的颜色是否小于l，如果是，则跳过当前循环的剩余部分，继续下一次循环。*/
        /*调用intersect_simd4x函数，计算与顶点u相邻的顶点与子图顶点数组的交集，并将结果存储在sub_v[l - 1]中，同时更新sub_v_size[l - 1]的值。*/
        sub_v_size[l - 1] = intersect_simd4x(sub_v[l], sub_v_size[l], DAG_adj[u], DAG_deg[0][u], sub_v[l - 1]);

        if (l == 2) {/*判断变量l是否等于2。*/
            for (j = 0; j < sub_v_size[l - 1]; j++) { /*如果l等于2，则执行这个循环，遍历交集数组sub_v[l - 1]的所有元素。*/
                (*cliques)++; /* 对cliques指针指向的值进行自增操作。*/
            }
        }

        else { /*如果变量l不等于2，则执行这个else语句块。*/
            if (sub_v_size[l - 1] >= l - 1) { /*判断交集数组sub_v[l - 1]的大小是否大于等于l - 1。*/
                EBBkC_plus(l - 1, cliques); /* 如果条件满足，则调用EBBkC_plus函数，参数为l - 1和cliques。*/
            }
        }
    }
}


void EBBkC_Graph_t::EBBkC_plus_plus(int l, unsigned long long *cliques) {
    int c, i, j, k, p, e, e_, u, v, w, s, t, end, dist; /*声明了一系列整数变量，这些变量用于后续的循环和条件判断。*/
    /*判断子图的顶点数组大小sub_v_size[l]是否小于l，或者子图的边数组大小sub_e_size[l]是否小于l * (l - 1) / 2。*/
    if (sub_v_size[l] < l || sub_e_size[l] < l * (l - 1) / 2) return; /*如果任一条件满足，则函数直接返回。*/

    if (l == K) { /*判断变量l是否等于类成员变量K。*/
        if (K == 3) { /*如果K等于3，则执行下面的代码块。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*循环遍历子图的边数组sub_e[l]的所有元素。*/
                e = sub_e[l][i]; /*获取子图的边数组的第i个元素，并将其存储在变量e中。*/
                for (j = 0; j < T_size[e]; j++) { /*对每条边e，循环遍历与该边相关的数组T[e]的所有元素。*/
                    (*cliques)++; /*对指针cliques指向的值进行自增操作。*/
                }
            }
        }
        else if (K == 4) { /*如果上面的条件不满足且K等于4，则执行下面的代码块。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*循环遍历子图的边数组sub_e[l]的所有元素。*/
                e = sub_e[l][i]; /*获取子图的边数组的第i个元素，并将其存储在变量e中。*/
                for (j = 0; j < C_size[e]; j++) { /* 对每条边e，循环遍历与该边相关的数组C[e]的所有元素。*/
                    (*cliques)++; /*对指针cliques指向的值进行自增操作。*/
                }
            }
        }
        else {/*之前的代码段处理了K等于3和4的特殊情况，这个else语句块则处理其他情况，即K不等于3或4的情况。*/
            for (i = 0; i < sub_e_size[l]; i++) { /*循环遍历子图的边数组sub_e[l]的所有元素。*/
                e = sub_e[l][i]; /*获取子图的边数组的第i个元素，并将其存储在变量e中。*/
                /*判断与边e相关的两个数组T[e]和C[e]的大小是否满足条件，如果不满足则跳过当前循环的剩余部分。*/
                if (T_size[e] < l - 2 || C_size[e] < (l - 2) * (l - 3) / 2) continue;
                
                for (j = 0; j < T_size[e]; j++) { /*对每条边e，循环遍历与该边相关的数组T[e]的所有元素。*/
                    u = T[e][j]; /*获取数组T[e]的第j个元素，并将其存储在变量u中。*/
                    col[u] = 0; /*将顶点u的颜色设置为0。*/
                    DAG_deg[l - 2][u] = 0; /* 将顶点u在DAG（有向无环图）中的度设置为0。*/
                    G_deg[l - 2][u] = 0; /*将顶点u在图G中的度设置为0。*/
                }

                for (j = 0; j < C_size[e]; j++) { /*对每条边e，循环遍历与该边相关的数组C[e]的所有元素。*/
                    e_ = C[e][j]; /*获取数组C[e]的第j个元素，并将其存储在变量e_中。*/
                    s = edges[e_].s; /*获取边e_的起点，并将其存储在变量s中。*/
                    t = edges[e_].t; /*获取边e_的终点，并将其存储在变量t中。*/
                    G_adj[s][G_deg[l - 2][s]++] = t; /*在图G的邻接矩阵中，将顶点s和顶点t连接起来，并更新顶点s的度。*/
                    G_adj[t][G_deg[l - 2][t]++] = s; /*在图G的邻接矩阵中，将顶点t和顶点s连接起来，并更新顶点t的度。*/
                }

                auto *list = new KeyVal_t [truss_num + 1]; /*: 动态分配了一个KeyVal_t类型的数组list，大小为truss_num + 1。*/
                for (j = 0; j < T_size[e]; j++) { /*循环遍历与边e相关的数组T[e]的所有元素。*/
                    u = T[e][j]; /*获取数组T[e]的第j个元素，并将其存储在变量u中。*/
                    list[j].key = u; /*将顶点u存储在list数组的第j个元素的键（key）字段中。*/
                    list[j].val = G_deg[l - 2][u]; /*获取顶点u在图G中的度，并将其存储在list数组的第j个元素的值（val）字段中。*/
                }
                sort(list, list + T_size[e]); /*对list数组进行排序，排序的依据是键值对中的键（key）。*/

                for (j = 0; j < T_size[e]; j++) { /*再次循环遍历与边e相关的数组T[e]的所有元素，但这次是通过已排序的list数组进行。*/
                    u = list[j].key; /*获取已排序的list数组的第j个元素的键（key）字段，并将其存储在变量u中。*/
                    for (k = 0; k < G_deg[l - 2][u]; k++) { /*循环遍历与顶点u相邻的所有顶点。*/
                        v = G_adj[u][k]; /*获取顶点u的第k个相邻顶点，并将其存储在变量v中。*/
                        used[K][col[v]] = true; /*将二维数组used的对应位置设置为true，表示该顶点已被使用或访问过。*/
                    }
                    for (c = 1; used[K][c]; c++) ; /*查找第一个未被使用的颜色。*/
                    col[u] = c; /*将顶点u的颜色设置为找到的第一个未被使用的颜色。*/
                    for (k = 0; k < G_deg[l - 2][u]; k++) { /*再次循环遍历与顶点u相邻的所有顶点。*/
                        v = G_adj[u][k]; /*获取顶点u的第k个相邻顶点，并将其存储在变量v中。*/
                        used[K][col[v]] = false; /*将二维数组used的对应位置设置为false，表示该颜色已被使用，不能被其他顶点使用。*/
                    }
                }
                delete [] list; /*删除动态分配的数组list，释放内存空间。*/

                sub_v_size[l - 2] = 0; /*将子图的顶点数组sub_v[l - 2]的大小重置为0。*/
                dist = 0; /* 初始化变量dist为0，用于记录具有不同颜色的顶点数量。*/

                for (j = 0; j < T_size[e]; j++) { /*循环遍历与边e相关的数组T[e]的所有元素。*/
                    u = T[e][j]; /*获取数组T[e]的第j个元素，并将其存储在变量u中。*/
                    sub_v[l - 2][sub_v_size[l - 2]++] = u; /* 将顶点u添加到子图的顶点数组sub_v[l - 2]中，并更新数组大小。*/
                    if (!used[K][col[u]]) { /*检查顶点u的颜色是否已被使用过。*/
                        used[K][col[u]] = true; /*如果颜色未被使用，将其标记为已使用。*/
                        dist++; /*增加变量dist的值，表示找到一个具有不同颜色的顶点。*/
                    }
                }

                if (dist >= l - 2) { /*检查具有不同颜色的顶点数量是否大于等于l - 2。*/
                    sub_e_size[l - 2] = 0; /*将子图的边数组sub_e[l - 2]的大小重置为0。*/
                    for (j = 0; j < C_size[e]; j++) { /*循环遍历与边e相关的数组C[e]的所有元素。*/
                        e_ = C[e][j]; /*获取数组C[e]的第j个元素，并将其存储在变量e_中。*/
                        sub_e[l - 2][sub_e_size[l - 2]++] = e_; /*将边e_添加到子图的边数组sub_e[l - 2]中，并更新数组大小。*/
                        s = edges[e_].s; /*获取边e_的起点，并将其存储在变量s中。 */
                        t = edges[e_].t; /*获取边e_的终点，并将其存储在变量t中。*/
                        edges[e_].s = (col[s] > col[t]) ? s : t; /*它检查s的颜色是否大于t的颜色。如果col[s]大于col[t]，则边e_的起点保持不变为s；否则，将边e_的起点设置为t。*/
                        edges[e_].t = (col[s] > col[t]) ? t : s; /*与上一句相反，它检查s的颜色是否小于或等于t的颜色。如果col[s]大于col[t]，则边e_的终点设置为t；否则，边e_的终点设置为s*/
                        s = edges[e_].s; /*重新获取边e_更新后的起点，并将其存储在变量s中。*/
                        t = edges[e_].t; /*重新获取边e_更新后的终点，并将其存储在变量t中。*/

                        DAG_adj[s][DAG_deg[l - 2][s]++] = t; /*在有向无环图（DAG）的邻接矩阵中添加从顶点s到顶点t的边，并更新顶点s的度。*/
                    }

                    EBBkC_plus_plus(l - 2, cliques); /*这行代码是一个递归调用，调用了函数EBBkC_plus_plus并将l - 2和cliques作为参数传递给它。*/
                }

                for (j = 0; j < T_size[e]; j++) { /*循环，用于重置used数组中的某些值。*/
                    u = T[e][j]; /*循环遍历与边e相关的数组T[e]的所有元素，*/
                    used[K][col[u]] = false; /*并将每个元素对应的颜色在used数组中的值设置为false，表示该颜色已被释放或未被使用。*/
                }
            }
        }

        return;
    }

    if (l == 2) { /*如果 l 等于 2，则执行以下的代码块  */
        for (i = 0; i < sub_v_size[l]; i++) { /*循环遍历子图的所有顶点。*/
            u = sub_v[l][i]; /*从sub_v数组中获取对应的顶点u  */

            for (j = 0; j < DAG_deg[l][u]; j++) { /*对于顶点u，开始一个循环，从0开始，直到小于DAG_deg[l][u]的值*/
                (*cliques)++; /*cliques指针指向的值增加1，表示找到一个团或子结构  */
            }
        }

        return;
    }

    if (l == 3) { /*当 l 等于 3 时，执行以下的代码块  */

        for (i = 0; i < sub_v_size[l]; i++) { /*遍历子图的所有顶点 */
            u = sub_v[l][i]; /*获取当前子图的顶点 u  */

            if (col[u] < l) continue; /*如果顶点 u 的颜色小于 l，则跳过此次循环*/

            for (j = 0; j < DAG_deg[l][u]; j++) { /*遍历 u 在 DAG 中的邻居。*/
                v = DAG_adj[u][j]; /*获取 u 的邻居 v */
                lab[v] = l - 1; /*将 v 的标签设置为 l-1  */
            }

            for (j = 0; j < DAG_deg[l][u]; j++) { /*再次遍历 u 在 DAG 中的邻居 */
                v = DAG_adj[u][j]; /*获取 u 的邻居 v 。*/

                if (col[v] < l - 1) continue; /*如果顶点 v 的颜色小于 l-1，则跳过此次循环。*/

                for (k = 0; k < DAG_deg[l][v]; k++) { /*遍历 v 在 DAG 中的邻居。*/
                    w = DAG_adj[v][k]; /*获取 v 的邻居 w。*/
                    if (lab[w] == l - 1) (*cliques)++; /*如果 w 的标签为 l-1，则团的数量增加1。*/
                }
            }

            for (j = 0; j < DAG_deg[l][u]; j++) { /*最后遍历一次 u 在 DAG 中的邻居。*/
                v = DAG_adj[u][j]; /*获取 u 的邻居 v*/
                lab[v] = l; /*将 v 的标签重置为 l  */
            }
        }

        return;
    }

    if (can_terminate(l, cliques)) { /*如果可以终止条件满足，则执行以下代码  */
        return; /*直接返回，结束当前函数。*/
    }

    for (i = 0; i < sub_v_size[l]; i++) { /*遍历子图的所有顶点。*/
        u = sub_v[l][i]; /*获取当前子图的顶点 u  */

        if (col[u] < l) continue; /*如果顶点 u 的颜色小于 l，则跳过此次循环  */

        sub_v_size[l - 1] = 0; /*将子图 l-1 的大小重置为0 */
        dist = 0; /*初始化距离为0  */

        for (j = 0; j < DAG_deg[l][u]; j++) { /*遍历 u 在 DAG 中的邻居  */
            v = DAG_adj[u][j]; /*获取 u 的邻居 v  */
            lab[v] = l - 1; /*将 v 的标签设置为 l-1  */
            sub_v[l - 1][sub_v_size[l - 1]++] = v; /*将 v 加入到子图 l-1 中，并更新子图 l-1 的大小*/
            DAG_deg[l - 1][v] = 0; /*将 v 在 DAG 中 l-1 层的度重置为0*/
            G_deg[l - 1][v] = 0; /*将 v 在 G 中 l-1 层的度重置为0  */

            if (!used[l][col[v]]) { /*如果 v 的颜色在 used 数组中的对应位置为 false  */
                used[l][col[v]] = true; /*则将 used 数组中 v 的颜色对应位置设置为 true  */
                dist++; /*距离增加1  */
            }
        }

        if (dist >= l - 1) { /*如果距离大于或等于 l-1，则执行以下代码块  */

            sub_e_size[l - 1] = 0; /*将子图 l-1 的边的大小重置为0  */
            for (j = 0; j < sub_v_size[l - 1]; j++) { /*遍历子图 l-1 的所有顶点  */
                v = sub_v[l - 1][j]; /*获取当前子图 l-1 的顶点 v  */

                end = DAG_deg[l][v]; /*获取 v 在 DAG 中 l 层的度，并将其赋值给 end */
                for (k = 0; k < end; k++) { /*遍历 v 在 DAG 中的邻居  */
                    w = DAG_adj[v][k]; /*获取 v 的邻居 w  */
                    if (lab[w] == l - 1) { /*如果 w 的标签等于 l-1 。*/
                        DAG_deg[l - 1][v]++; /*将 v 在 DAG 中 l-1 层的度增加1  */
                        sub_e_size[l - 1]++; /*子图 l-1 的边的大小增加1 */

                        // just for early-termination
                        G_deg[l - 1][v]++; /*将 v 在 G 中 l-1 层的度增加1（仅用于早期终止）*/
                        G_deg[l - 1][w]++; /*将 w 在 G 中 l-1 层的度增加1（仅用于早期终止）*/

                    } else {
                        DAG_adj[v][k--] = DAG_adj[v][--end]; /*将 DAG_adj 数组中 v 的邻居 w 替换为其后面的邻居  */
                        DAG_adj[v][end] = w; /*将 w 放在 DAG_adj 数组中 v 的邻居的末尾  */
                    }
                }
            }

            EBBkC_plus_plus(l - 1, cliques); /*调用 EBBkC_plus_plus 函数，传入参数 l-1 和 cliques  */
        }

        for (j = 0; j < sub_v_size[l - 1]; j++) { /*遍历子图 l-1 的所有顶点  */
            v = sub_v[l - 1][j]; /*获取子图 l-1 的当前顶点 v */
            lab[v] = l; /*将 v 的标签设置为 l  */
            used[l][col[v]] = false; /*将 used 数组中，索引为 l、col[v] 的元素设置为 false*/
        }
    }
}


void EBBkC_Comb_list(int *list, int list_size, int start, int picked, int k, unsigned long long *cliques) { 
    /*用于计算从给定列表中选取k个元素的组合数。这是一种基于递归的方法。*/
    if (picked == k) { /*检查是否已选取的元素数量等于需要的元素数量。*/
        (*cliques)++; /*如果上述条件为真，则增加cliques所指向的值。这实际上是在计数符合条件的组合数。*/
        return;
    }

    for (int i = start; i < list_size; i++) { /*从start开始遍历列表，直到list_size。*/
        EBBkC_Comb_list(list, list_size, i + 1, picked + 1, k, cliques); /*这是函数的递归调用。每次循环时，它都会调用自身，但是会更新一些参数：
        start被设置为i + 1，意味着下一次递归调用将从当前索引的下一个位置开始。picked被增加1，表示已经选取了一个额外的元素。*/
    }
}

void EBBkC_Graph_t::list_in_plex(int start, int p, int q, unsigned long long *cliques) {
    if (F_size < q) return; /*如果F_size小于q，则直接返回，不执行后续操作 */

    if (p == 0) {
        if (q > F_size - q) 
            EBBkC_Comb_list(F, F_size, 0, 0, F_size - q, cliques); /*调用EBBkC_Comb_list函数*/
        else
            EBBkC_Comb_list(F, F_size, 0, 0, q, cliques); /*否则调用EBBkC_Comb_list函数，参数为F, F_size, 0, 0, q和cliques  */
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