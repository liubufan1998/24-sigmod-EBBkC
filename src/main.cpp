#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <set>
#include <algorithm>
#include <string>
#include <omp.h>
#include <cassert>
#include "set_operation.h"
#include "def.h"
#include "edge_oriented.h"


using namespace std;
int K, L;
unsigned long long N = 0;

/*这段代码是一个用于解决k-clique问题的程序，主要包含两个部分：预处理（pre-process）和EBBkC+ET算法执行（EBBkC+ET）。*/
int main(int argc, char** argv) {
    double runtime;
    string act(argv[1]);

    if (act == "p") { // pre-process。这部分的功能主要是对输入的边列表文件（.edges或.mtx格式）进行清理，然后进行网络的k-clique本征值分解，将结果存储在索引文件中。
        printf("Pre-process %s\n", argv[2]); /*打印出正在处理的文件名（argv[2]）,argv[2]是命令行参数数组中的第三个元素，它应该是输入文件的名称。*/
        string src_filename(argv[2]); /*创建一个字符串对象src_filename，并将其初始化为argv[2]的值，也就是输入文件的名称*/
        string suffix = src_filename.substr(src_filename.find_last_of('.')); /*从输入文件名中提取出后缀名。*/
        if (suffix != ".edges" && suffix != ".mtx") exit(0); /*判断输入的文件后缀名是否为".edges"或".mtx"，如果不是则退出程序。*/
        string prefix = src_filename.substr(0, src_filename.find_last_of('.')); /*从输入文件名中提取出前缀名。*/
        string clean_filename = prefix + ".clean"; /*构造一个新的字符串clean_filename，它是前缀名加上".clean"。这可能用于存储清理后的版本的文件。*/
        clean_edges(src_filename.c_str(), clean_filename.c_str()); /*调用一个名为clean_edges的函数,用于清理边列表文件*/
        string index_filename = prefix + ".index"; /*构造一个新的字符串index_filename，它是前缀名加上".index"。这可能用于存储索引文件*/
        runtime = EBBkC_t::truss_order(clean_filename.c_str(), index_filename.c_str()); /*调用一个名为EBBkC_t::truss_order的函数*/
        printf("Pre-processed in %.2lf ms\n\n", runtime); /*告知用户预处理所花费的时间。*/
    }

    else if (act == "e") { // EBBkC+ET。读取预处理生成的索引文件，并使用EBBkC+ET算法来寻找网络中的k-clique
        string src_filename(argv[2]); /*初始化一个字符串对象src_filename，将其值设置为argv[2]。argv[2]应该是命令行参数数组中的第三个元素，它代表输入文件的名称。*/
        string suffix = src_filename.substr(src_filename.find_last_of('.')); /*从输入文件名中提取后缀名。*/
        if (suffix != ".index") exit(0); /*检查提取出的后缀名是否为".index"。如果不是，则程序退出。*/

        K = atoi(argv[3]); /*将命令行参数数组中的第四个元素转换为整数，并将其赋值给变量K。K是k-clique的k值。*/
        L = atoi(argv[4]); /*将命令行参数数组中的第五个元素转换为整数，并将其赋值给变量L.*/

        runtime = EBBkC_t::list_k_clique(argv[2]); /*调用一个名为list_k_clique的函数，并将运行时间存储在变量runtime中.*/
        printf("Number of %u-cliques: %llu\n", K, N); /*打印出找到的k-clique的数量。*/
        printf("EBBkC+ET (t = %d) runtime %.2lf ms\n\n", L, runtime); /*打印出EBBkC+ET算法的运行时间。*/
    }

    else { /*否则输入错误，提示用户错误信息。*/
        printf("Wrong usage.\n");
        exit(0);
    }

    return 0;
}
