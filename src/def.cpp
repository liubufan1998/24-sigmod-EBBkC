#include <cassert>
#include "def.h"
#include <set>

size_t h(const Edge_t& e) { /*定义了一个名为h的函数，它接受一个Edge_t类型的常量引用作为参数，并返回一个size_t类型的值。该函数将用于计算边的哈希值。*/
    int s_ = e.s < e.t ? e.s : e.t; /*使用三元运算符来判断边e的两个节点s和t的大小，并将较小的值赋给变量s_.*/
    int t_ = e.s < e.t ? e.t : e.s; /*使用三元运算符来判断边e的两个节点s和t的大小，并将较大的值赋给变量t_。*/
    size_t hash = 1; /*该变量将用于存储计算出的哈希值。*/
    std::hash<int> H; /*该对象将用于计算整数的哈希值。*/

    hash ^= H(s_) + 0x9e3779b9 + (hash << 6) + (hash >> 2); /*首先对s_调用哈希函数，然后与一些常数和hash的位移版本进行异或操作。*/
    hash ^= H(t_) + 0x9e3779b9 + (hash << 6) + (hash >> 2); /*这行与上一行类似，但是它对t_进行哈希计算，并再次更新hash的值。*/
    return hash; /*返回计算出的哈希值。*/
}

Edge_t::Edge_t() = default; /*类的默认构造函数的定义，它使用默认实现。*/

Edge_t::Edge_t(int s, int t, bool directed) { /*定义了一个构造函数，它接受两个整数和一个布尔值作为参数。这两个整数表示边的两个节点，布尔值表示边是否有向。*/
    if (directed) { /*如果边是有向的。*/
        this->s = s; /*将参数s和t分别赋给类的成员变量s和t*/
        this->t = t;
    }
    else { /*如果边是无向的，这行代码将使用三元运算符来确保节点按照升序存储，即较小的节点存储在s中，较大的节点存储在t中。*/
        this->s = s < t ? s : t;
        this->t = s < t ? t : s;
    }
}

bool Edge_t::operator<(const Edge_t &e) const { /*定义了一个小于运算符的重载，用于比较两个Edge_t对象。*/
    return (this->s < e.s) || (this->s == e.s && this->t < e.t); /*它首先比较成员变量s，如果相等则比较成员变量t。这使得对象可以按照字典序进行排序。*/
}

bool Edge_t::operator==(const Edge_t &e) const { /*定义了一个等于运算符的重载，用于判断两个Edge_t对象是否相等。*/
    return this->s == e.s && this->t == e.t; /*它通过比较成员变量s和t来判断对象是否相等。*/
}

size_t Edge_t::Hash_Edge_t::operator()(const Edge_t &e) const { /*定义了一个名为Hash_Edge_t的嵌套类的函数调用运算符的重载*/
    return h(e); /*该类用于计算边的哈希值，并通过调用前面定义的函数h来实现这一点。*/
}

HashMap_t::HashMap_t() { /*定义了一个名为HashMap_t的哈希表类，用于存储和查找关于K-Clique问题的边的信息。定义了HashMap_t类的构造函数，用于初始化哈希表。*/
    table = new HashMap_t::HashItem_t* [TAB_SIZE]; /*这个数组将用作哈希表的存储结构，其中每个元素指向一个链表，用于解决哈希冲突。*/
    table_size = new int [TAB_SIZE](); /*这个数组用于记录每个哈希桶中元素的数量。*/
}

HashMap_t::~HashMap_t() { /*定义了HashMap_t类的析构函数，用于释放哈希表占用的内存。*/
    delete [] table;
    delete [] table_size;
}

int HashMap_t::exist(const Edge_t &key) { /*用于检查给定的边是否已经存在于哈希表中。*/
    size_t pos = h(key) % TAB_SIZE; /*计算了给定边的哈希值，并通过取模运算将其映射到哈希表的一个特定位置。这里使用了前面定义的哈希函数h。*/
    int i, val; /*变量i用于循环计数, 变量val用于存储找到的边的值。*/
    Edge_t e; /*用于临时存储哈希表中的边。*/

    if (table_size[pos] > 0) { /*如果桶的大小大于0，表示有元素需要比较。*/
        for (i = 0; i < table_size[pos]; i++) { /*遍历指定位置的哈希桶中的所有元素。循环的次数由桶的大小决定。*/
            e = table[pos][i].key; /*从哈希桶中取出当前元素的键（即边），并将其存储在局部变量e中。*/
            val = table[pos][i].val; /*从哈希桶中取出当前元素的值，并将其存储在局部变量val中。这里的值可能是与边相关联的某种信息。*/
            /*检查当前元素是否与给定的边相等。由于边可能是无向的，所以需要同时比较两种可能的节点顺序。*/
            if ((e.s == key.s && e.t == key.t) || (e.s == key.t && e.t == key.s)) return val; /*如果找到了匹配的边，函数将返回其关联的值。*/
        }
    }

    return -1; /*如果没有找到匹配的边，则返回-1.*/
}

void HashMap_t::insert(const Edge_t &key, int val) { /*用于向哈希表中插入边。*/
    size_t pos = h(key) % TAB_SIZE; /*计算边key的哈希值，并通过取模运算确定其在哈希表中的位置. */

    if (table_size[pos] == 0) /*如果指定位置的哈希桶为空，*/
        table[pos] = new HashMap_t::HashItem_t [MAX_COLL]; /*则分配一个新的数组来存储哈希项。*/
    /*使用断言确保哈希桶的大小不会超过MAX_COLL。如果超过，则需要重新定义哈希函数或增大MAX_COLL的值。*/
    assert(table_size[pos] < MAX_COLL);     // Enlarge max_coll or redefine hash function h(e)

    table[pos][table_size[pos]].key = key; /*在哈希桶的末尾插入新的边。*/
    table[pos][table_size[pos]++].val = val; /*在哈希桶的末尾插入与边关联的值，并递增桶的大小。*/
}

void HashMap_t::remove(const Edge_t &key) { /*用于向哈希表中删除边, 该函数接受一个类型为Edge_t的常量引用参数key。*/
    size_t pos = h(key) % TAB_SIZE; /*计算边key的哈希值，并通过取模运算确定其在哈希表中的位置。*/
    int i; /*定义一个局部整数变量i用于循环计数*/
    Edge_t e; /*定义一个局部变量e用于临时存储哈希表中的边。*/

    if (table_size[pos] > 0) { /*检查指定位置的哈希桶是否为空。如果不为空，则执行删除操作。*/
        for (i = 0; i < table_size[pos]; i++) { /*遍历指定位置的哈希桶中的所有元素。*/
            e = table[pos][i].key; /*从哈希桶中取出当前元素的键（即边），并将其存储在局部变量e中。*/
            if ((e.s == key.s && e.t == key.t) || (e.s == key.t && e.t == key.s)) { /*检查当前元素是否与给定的边相等。由于边可能是无向的，所以需要同时比较两种可能的节点顺序。*/
                table[pos][i] = table[pos][--table_size[pos]]; /*如果找到了匹配的边，则用哈希桶中的最后一个元素覆盖当前元素，并递减桶的大小。这实际上删除了指定的边。*/
                if (table_size[pos] == 0) delete [] table[pos]; /*如果哈希桶变为空，则释放其占用的内存。*/
                break; /*跳出循环，因为已经找到了并删除了指定的边。*/
            }
        }
        return; /*结束函数的执行，返回调用者。如果边不存在于哈希表中，则直接返回。*/
    }

    exit(-1); /*如果代码执行到这里，表示给定的边不存在于哈希表中。这一行将终止程序的执行，并返回一个错误代码-1。*/
}

/*下面的这段代码定义了两个类：KeyVal_t 和 Heap_t，其中 KeyVal_t 是一个用于存储键值对的简单结构，
而 Heap_t 则是一个基于 KeyVal_t 的堆结构，用于维护键值对的最大堆。*/
KeyVal_t::KeyVal_t() = default; /* KeyVal_t 类的默认构造函数的定义，使用 = default; 指示编译器自动生成默认的实现。*/

KeyVal_t::KeyVal_t(int key, int val) { /*KeyVal_t 类的带参数构造函数的定义，接受两个整数参数 key 和 val。*/
    this->key = key; /*将传入的 key 参数值赋给成员变量 key。*/
    this->val = val;
}

bool KeyVal_t::operator<(const KeyVal_t &kv) const { /*定义了一个小于运算符的重载函数，用于比较两个 KeyVal_t 对象。*/
    return this->val > kv.val || (this->val == kv.val && this->key < kv.key); /*比较规则是先按 val 的值降序排列，如果 val 相等，则按 key 的值升序排列。*/
}

Heap_t::Heap_t() = default; /*Heap_t 类的默认构造函数的定义，使用 = default; 指示编译器自动生成默认的实现。*/
Heap_t::~Heap_t() = default; /*Heap_t 类的析构函数的定义，使用 = default; 指示编译器自动生成默认的实现。*/

void Heap_t::swap(unsigned i, unsigned j) { /*定义了一个名为 swap 的成员函数，用于交换堆中两个位置的元素。*/
    KeyVal_t kv_tmp = this->kv_list[i]; /*创建一个临时的 KeyVal_t 对象 kv_tmp 并将其初始化为堆列表中位置 i 的元素。*/
    int pt_tmp = this->pt[kv_tmp.key]; /*创建一个临时的整数变量 pt_tmp 并将其初始化为与 kv_tmp.key 相关联的值*/
    this->pt[this->kv_list[i].key] = this->pt[this->kv_list[j].key]; /*将位置 j 的元素的 pt 值赋给位置 i 的元素的 pt 值。*/
    this->kv_list[i] = this->kv_list[j]; /*将位置 j 的元素赋给位置 i 的元素。*/
    this->pt[this->kv_list[j].key] = pt_tmp; /*将临时变量 pt_tmp 的值赋给位置 j 的元素的 pt 值。*/
    this->kv_list[j] = kv_tmp; /*将临时变量 kv_tmp 的值赋给位置 j 的元素。*/
    /*这个 swap 函数不仅交换了堆中的两个元素，还更新了与这两个元素相关的 pt 数组的值。这样可以确保堆和 pt 数组之间的一致性。*/
}

/*下面这段代码是Heap_t类的两个成员函数：bubble_up和bubble_down。这两个函数是用于维护最大堆性质的。*/
void Heap_t::bubble_up(unsigned int i) { /*这个函数是为了确保新加入的元素在堆中正确地找到自己的位置，从而保证堆的性质。*/
    unsigned j = (i - 1) >> 1; /*计算元素i的父节点的索引，并将其存储在j中。这里使用了位操作(i - 1) >> 1，这是一种高效的计算父节点索引的方法。*/
    while (i > 0) { /*开始一个循环，只要i不是根节点（即0），就继续循环。*/
        if (this->kv_list[j].val > this->kv_list[i].val) { /*检查元素i的值是否小于其父节点j的值。如果是，则交换这两个元素的位置。*/
            this->swap(i, j); /*交换元素i和元素j的位置。*/
            i = j; /*更新i的值为j，因为我们现在要考虑的是原来位置j的元素（现在已经在位置i上了）。*/
            j = (i - 1) >> 1; /*重新计算元素i的父节点的索引，并将其存储在j中。*/
        }
        else break; /*如果元素i的值不小于其父节点j的值，那么堆的性质已经得到了维护，循环结束。*/
    }
}
/*个函数是为了确保被修改（例如减少键值）的元素在堆中正确地找到自己的位置，从而保证堆的性质。*/
void Heap_t::bubble_down(unsigned int i) { /*它接受一个无符号整数i作为参数，表示要下沉的元素的索引。*/
    unsigned l = (i << 1) + 1, r = l + 1, j; /*计算元素i的左子节点和右子节点的索引，并分别存储在l和r中。
    这里使用了位操作(i << 1) + 1和l + 1，这是一种高效的计算子节点索引的方法。*/
    while (l < this->n) {/*开始一个循环，只要左子节点存在（即其索引小于堆的大小），就继续循环。*/
        /*确定要与其父节点比较的子节点的索引。如果右子节点存在且其值小于左子节点的值，则选择右子节点，否则选择左子节点。*/
        j = ((r < this->n) && (this->kv_list[r].val < this->kv_list[l].val)) ? r : l;
        if (this->kv_list[i].val > this->kv_list[j].val) { /*检查元素i的值是否大于子节点j的值。如果是，则交换这两个元素的位置。*/
            this->swap(i, j); /*交换元素i和元素j的位置。*/
            i = j; /*更新i的值为j，因为我们现在要考虑的是原来位置j的元素（现在已经在位置i上了）。*/
            l = (i << 1) + 1; /*重新计算元素i的左子节点的索引，并将其存储在l中。*/
            r = l + 1; /*重新计算元素i的右子节点的索引，并将其存储在r中。*/
        }
        else break; /*如果元素i的值不大于其子节点的值，那么堆的性质已经得到了维护，循环结束。*/
    }
}
/*下面的这段代码是Heap_t类的几个成员函数，包括empty、insert、pop、min_element和update。这些函数分别用于检查堆是否为空、插入新元素、弹出最小元素、
获取最小元素和更新元素的值。*/
bool Heap_t::empty() { /*定义了一个名为empty的成员函数，用于检查堆是否为空。*/
    return this->n == 0; /*如果堆的大小（即成员变量n）为0，则返回true，表示堆为空；否则返回false。*/
}

void Heap_t::insert(KeyVal_t kv) { /*定义了一个名为insert的成员函数，接受一个KeyVal_t类型的参数kv，用于向堆中插入新元素。*/
    this->pt[kv.key] = this->n; /*在pt数组中记录新元素的键和它在堆中的索引。*/
    this->kv_list[this->n] = kv; /*将新元素添加到堆的末尾。*/
    this->bubble_up(this->n); /*通过调用bubble_up函数，确保新插入的元素在堆中正确地找到自己的位置，从而维护堆的性质。*/
    this->n++; /*增加堆的大小。*/
}

KeyVal_t Heap_t::pop() { /*定义了一个名为pop的成员函数，用于从堆中弹出并返回最小元素。*/
    assert(!this->empty()); /*使用断言确保堆不为空，如果为空则触发错误。*/
    KeyVal_t min = this->kv_list[0]; /*获取堆顶元素（即最小元素）并存储在变量min中。*/
    this->pt[min.key] = -1; /*在pt数组中将弹出的元素的键对应的索引设置为-1，表示该键不再在堆中。*/
    this->kv_list[0] = this->kv_list[--(this->n)]; /*将堆的最后一个元素移到堆顶位置，并减少堆的大小。*/
    this->pt[this->kv_list[0].key] = 0; /*更新pt数组中新的堆顶元素的索引。*/
    this->bubble_down(0); /*通过调用bubble_down函数，确保替换后的堆顶元素在堆中正确地找到自己的位置，从而维护堆的性质。*/
    return min; /*返回弹出的最小元素。*/
}

KeyVal_t Heap_t::min_element() { /*定义了一个名为min_element的成员函数，用于返回堆中的最小元素（即堆顶元素）。*/
    assert(!this->empty()); /*使用断言确保堆不为空，如果为空则触发错误。*/
    return this->kv_list[0]; /*返回堆顶元素（即最小元素）。*/
}

void Heap_t::update(unsigned int key) { /*定义了一个名为update的成员函数，接受一个无符号整数参数key，用于更新堆中指定键的元素的值。*/
    int i = this->pt[key]; /*通过键在pt数组中查找元素的索引。*/
    if (i != -1) { /*检查找到的索引是否有效（即不为-1）。*/
        ((this->kv_list[i]).val)--; /*减少找到的元素的值。这里使用了复杂的语法来减少对成员变量val的引用级别。*/
        this->bubble_up(i); /*通过调用bubble_up函数，确保更新后的元素在堆中正确地找到自己的位置，从而维护堆的性质。*/
    }
}

/*下面这段代码是Heap_t类的两个成员函数，make_heap和release_heap。这些函数分别用于从一个给定的值列表创建一个堆，并释放堆所使用的内存*/
void Heap_t::make_heap(const int *v_list, unsigned int v_size) { /*定义了一个名为make_heap的成员函数，它接受一个整数指针v_list和一个无符号整数
v_size作为参数。这个函数用于从一个给定的整数列表中创建一个堆。*/
    unsigned i; /*声明一个无符号整数变量i，用于循环控制。*/
    KeyVal_t kv; /*声明一个KeyVal_t类型的变量kv，用于临时存储要插入堆的键值对。*/

    this->n = 0; /*将堆的大小（即成员变量n）设置为0，表示堆为空。*/
    this->pt = new int [v_size]; /*这个数组用于存储每个键在堆中的索引。*/
    for (i = 0; i < v_size; i++) this->pt[i] = -1; /*初始化pt数组，将所有元素的值设置为-1，表示目前没有任何键在堆中。*/
    this->kv_list = new KeyVal_t [v_size]; /*这个数组用于存储堆中的键值对。*/

    for (i = 0; i < v_size; i++) { /*开始一个循环，遍历给定的整数列表。*/
        kv.key = i; /*设置键值对的键为当前的索引。*/
        kv.val = v_list[i]; /*设置键值对的值为当前索引对应的值从给定的整数列表中。*/
        this->insert(kv); /*调用之前定义的insert函数，将当前的键值对插入到堆中。*/
    }
}

void Heap_t::release_heap() { /*定义了一个名为release_heap的成员函数，用于释放堆所使用的内存。*/
    delete [] this->kv_list; /*释放成员变量kv_list指向的动态分配的数组内存。*/
    delete [] this->pt; /*释放成员变量pt指向的动态分配的数组内存。*/
}

void clean_edges(const char *r_file_name, const char *w_file_name) {/*接受两个参数：r_file_name（输入文件的名称）和w_file_name（输出文件的名称）。*/
    set<Edge_t> E; /*这个集合用来存储读取的边，以消除重复的边。*/
    int s, t; /*用来存储每条边的两个顶点。*/
    FILE *f_r = fopen(r_file_name, "r"); /*打开一个文件，文件名由参数r_file_name指定，以只读模式（"r"）。*/
    FILE *f_w = fopen(w_file_name, "w"); /*打开另一个文件，文件名由参数w_file_name指定，以写入模式（"w"）*/

    while (fscanf(f_r, "%u %u%*[^\n]%*c", &s, &t) == 2) {/*从文件f_r中读取边。每行的格式应该是两个整数（用空格隔开），后面可能跟着一些非换行字符，
    然后是一个换行符。这个循环会持续读取文件的下一行，直到文件的末尾。*/
        if (s == t) {/*检查读取的边是否是一个自环（即两个顶点相同）*/
            printf("%u Self loop.\n", s); /*如果是自环，输出一个消息并跳过当前循环。*/
            continue;
        }
        if (E.find(Edge_t(s, t, false)) != E.end()) { /*检查读取的边是否已经在集合E中。如果找到相同的边，跳到下一个循环。*/
            printf("(%u, %u) duplicates.\n", s, t); /*如果找到重复的边，输出一个消息。*/
            continue;
        }

        E.insert(Edge_t(s, t, false)); /*将新的边（从文件读取的）插入集合E.使用三个参数构造一个Edge_t对象：两个顶点和一个标志位表示这个边没有被处理过。*/
        fprintf(f_w, "%u %u\n", s, t); /*将处理过的边写入文件f_w。每个边由两个顶点组成，后面跟着一个换行符。*/
    }

    fclose(f_r); /*关闭输入文件。*/
    fclose(f_w); /*关闭输出文件。*/
}

void GetCurTime(struct rusage* curTime) { /*获取当前的时间使用情况，并将其保存在传入的 curTime 指针指向的 rusage 结构体中。*/
    if (getrusage(RUSAGE_THREAD, curTime) != 0) { /*调用 getrusage 函数来获取当前线程的资源使用情况。如果该函数调用失败，它会返回一个非零值。*/
        fprintf(stderr, "The running time info couldn't be collected successfully.\n"); /*如果 getrusage 函数调用失败，这行代码会将错误信息打印到标准错误输出。*/
        exit(0); /*在打印错误信息后，程序会立即退出，返回值为0。*/
    }
}

/*该函数的目的是计算两个时间点之间的时间差，单位为毫秒。这两个时间点通过两个 rusage 结构体指针 start 和 end 传入。*/
double GetTime(struct rusage* start, struct rusage* end) { /*计算两个时间点之间的秒数差，并将结果转换为毫秒。这里，ru_utime.tv_sec 表示用户模式下的CPU时间.
计算两个时间点之间的微秒数差，并将结果转换为毫秒。这里，ru_utime.tv_usec 表示用户模式下的CPU时间（微秒）。这两部分的结果相加，得到两个时间点之间的总时间差，
单位为毫秒。*/
    return ((float)(end->ru_utime.tv_sec - start->ru_utime.tv_sec)) * 1e3 +
           ((float)(end->ru_utime.tv_usec - start->ru_utime.tv_usec)) * 1e-3;
}   // unit: ms

