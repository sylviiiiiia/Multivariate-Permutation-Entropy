#include<iostream>
#include<vector>
#include<algorithm>
#include<fstream>
#include<iomanip>
#include<cmath>
using namespace std;
#define int long long
const int N1 = 10, N2 = 10001;//时间序列数据,最大支持10通道，10001样本

//按second(数据)升序对pair排序
bool cmp(pair<int, double>& a, pair<int, double>& b) {
    return a.second < b.second;
}

//factorial计算阶乘
int factorial(int x) {
    int res = 1;
    for (int i = 2; i <= x; i++) {
        res *= i;
    }
    return res;
}

//argsort函数将a向量升序排序后的索引向量赋给t向量
void argsort(vector<double>& a, vector<int>& t, int len) {
    vector<pair<int, double> > nd(len);
    for (int i = 0; i < len; i++) {
        nd[i].first = i;
        nd[i].second = a[i];
    }
    sort(nd.begin(), nd.end(), cmp);
    for (int i = 0; i < len; i++) {
        t[i] = nd[i].first;
    }
}

//check函数检查两向量是否相同
bool check(vector<int> t, vector<int> p, int len) {
    bool f = 1;
    for (int i = 0; i < len; i++) {
        if (t[i] != p[i]) {
            f = 0;
            break;
        }
    }
    return f;
}

/*
mpe函数计算多变量时间序列的多变量排列熵。
-pe_channel:单通道排列熵
-pe_cross:交叉排列熵
-m:嵌入维度
-d:时间滞后
-n:时间序列样本长度
-e:时间序列样本通道数
*/
void mpe(vector<vector<double> >& data, vector<double>& pe_channel, double& pe_cross, int m, int d, int n, int e) {
    int t = n - d * (m - 1);
    int fact = factorial(m);
    vector<vector<int> > permutation(fact, vector<int>(m));
    vector<int> rangeGenerator(m);//用于生成全排列的原始数组，即range(m)
    for (int i = 0; i < m; i++)rangeGenerator[i] = i;
    for (int i = 0; i < fact; i++) {
        for (int j = 0; j < m; j++)permutation[i][j] = rangeGenerator[j];
        next_permutation(rangeGenerator.begin(), rangeGenerator.end());
    }
    vector<vector<int> >count(e, vector<int>(fact, 0));
    vector<vector<double> > possibility(e, vector<double>(fact));
    for (int s = 0; s < e; s++) {//第s个通道
        for (int k = 0; k < t; k++) {//处理第t个向量
            //对data按照k:k+d*m:d的方式切片,存入tmpdata
            vector<double> tmpdata(m);
            int j = 0;
            for (int i = k; i < k + d * m; i += d) {
                tmpdata[j++] = data[s][i];
            }
            //sortedIndex是tmpdata向量升序排序后的索引向量
            vector<int> sortedIndex(m);
            argsort(tmpdata, sortedIndex, m);
            for (int i = 0; i < fact; i++) {//sortedIndex与permutation匹配
                if (check(sortedIndex, permutation[i], m)) {
                    count[s][i]++;
                    break;
                }
            }
        }
        //计算各通道概率
        for (int i = 0; i < fact; i++) {
            possibility[s][i] = 1.0 * count[s][i] / (t * e);
            if (possibility[s][i]) {//当概率不为0时
                pe_channel[s] -= possibility[s][i] * log2(possibility[s][i]);
            }
        }
    }
    //计算交叉通道概率
    vector<double> rp(fact);
    for (int i = 0; i < fact; i++) {
        for (int s = 0; s < e; s++) {
            rp[i] += possibility[s][i];
        }
    }
    for (int i = 0; i < fact; i++) {
        if (rp[i])pe_cross -= rp[i] * log2(rp[i]);
    }

    //输出中间结果
    cout << "m=" << m << "  d=" << d << endl;
    cout << "permutation=" << endl;
    for (int i = 0; i < fact; i++) {
        for (int j = 0; j < m; j++)cout << permutation[i][j] << " ";
        cout << endl;
    }
    cout << "count=" << endl;
    for (int i = 0; i < e; i++) {
        for (int j = 0; j < fact; j++)cout << count[i][j] << " ";
        cout << endl;
    }
    cout << "possibility=" << endl;
    for (int i = 0; i < e; i++) {
        for (int j = 0; j < fact; j++)cout << setw(8) << possibility[i][j] << " ";
        cout << endl;
    }
}
signed main() {
    vector<vector<double> >data(N1, vector<double>(N2));
    int e = 1, n = 1000;//e是通道数，n是样本长度
    string fileName = "signal01.txt";
    int m = 3, d = 2;
    int ifput = 0;
    cout << "使用默认数据（signal01.txt，e=1,n=1000,m=3,d=2)请输入0，自行输入数据请输入1：";
    cin >> ifput;
    if (ifput == 1) {
        cout << "请输入txt文件名：";
        cin >> fileName;
        cout << "请输入通道数e及样本长度n：";
        cin >> e >> n;
        cout << "请输入嵌入维度m及时间延迟d：";
        cin >> m >> d;
    }
    else if (ifput != 0) {
        cout << "invalid input" << endl;
        return 0;
    }

    ifstream infile;
    infile.open(fileName, ios::in);
    if (!infile.is_open()) {
        cout << "读取文件失败" << endl;
        return 0;
    }
    int x = 0, y = 0;
    while (infile >> data[x][y]) {
        y++;
        if (y >= n) {
            x++;
            y = 0;
        }
    }
    vector<double> pe_channel(e);
    double pe_cross = 0;
    mpe(data, pe_channel, pe_cross, m, d, n, e);
    for (int i = 0; i < e; i++) {
        cout << "第" << i << "个通道的排列熵=" << setprecision(17) << pe_channel[i] << endl;
    }
    cout << "多通道交叉排列熵=" << setprecision(17) << pe_cross << endl;

    return 0;
}
