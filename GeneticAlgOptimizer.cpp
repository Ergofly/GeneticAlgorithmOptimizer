//
// Created by ergofly on 2022/12/4.
//
#include<iostream>
#include<cmath>
#include<vector>
#include<ctime>
#include <random>
#include <functional>

#define PI 3.1415926535897923384

//随机数生成引擎，生成0到1之间的随机数
std::default_random_engine e(time(nullptr));
std::uniform_real_distribution<double> u(0, 1);

double f1(double x1, double x2)
{
    return pow(x1, 2) + pow(x2, 2);
}

double f2(double x1, double x2)
{
    return std::abs(x1) + std::abs(x2) + std::abs(x1) * std::abs(x2);
}

double f3(double x1, double x2)
{
    return pow(x1, 2) + pow(x1 + x2, 2);
}

double f4(double x1, double x2)
{
    return -x2 * sin(sqrt(std::abs(x1))) - x2 * sin(sqrt(std::abs(x2)));
}

double f5(double x1, double x2)
{
    return pow(x1, 2) - 10 * cos(2 * PI * x1) + 10 + pow(x2, 2) - 10 * cos(2 * PI * x2) + 10;
}

double f6(double x1, double x2)
{
    return 20 * (1 - exp(-0.2 * sqrt((pow(x1, 2) + pow(x2, 2)) / (float) 2))) + exp(1) -
           exp((cos(2 * PI * x1) + cos(2 * PI * x2)) / (float) 2);
}

//函数类，每个函数是一个类
class OptFunction
{
private:
    double X1_BOUND[2]{0, 0};
    double X2_BOUND[2]{0, 0};
    std::function<double(double, double)> f;
public:
    int x1Size; //X1的编码长度
    int x2Size; //X2的编码长度
    double x1Field;
    double x2Field;

    OptFunction()
    {
        this->f = nullptr;
        this->x1Size = 0;
        this->x2Size = 0;
        this->x1Field = 0;
        this->x2Field = 0;
    }

    OptFunction(std::function<double(double, double)> f, double x1Low, double x1High, double x2Low, double x2High,
                double eps)
    {
        this->X1_BOUND[0] = x1Low;
        this->X1_BOUND[1] = x1High;
        this->X2_BOUND[0] = x2Low;
        this->X2_BOUND[1] = x2High;
        this->f = f;
        this->x1Size = ceil(log2((this->X1_BOUND[1] - this->X1_BOUND[0]) / eps));
        this->x2Size = ceil(log2((this->X2_BOUND[1] - this->X2_BOUND[0]) / eps));
        this->x1Field = this->X1_BOUND[1] - this->X1_BOUND[0];
        this->x2Field = this->X2_BOUND[1] - this->X2_BOUND[0];
    }

    double getX1LowBound() { return X1_BOUND[0]; }

    double getX2LowBound() { return X2_BOUND[0]; }

    double evaluate(double x1, double x2)
    {
        return this->f(x1, x2);
    }
};


//遗传算法求解无约束函数优化问题主类
class GeneticAlgOptimizer
{
private:
    int N, E, Times;
    double Pc, Pm, epsilon;

    //解码，返回第i个个体的函数值
    double decode(std::vector<std::vector<int>> &individual, int i)
    {
        //解码值x1,x2
        double x1 = 0, x2 = 0;
        //解码x1
        int tempSum = 0;
        for (int j = optF.x1Size - 1, k = 0; j >= 0; j--, k++)
        {
            tempSum += (1 << k) * individual[i][j];
        }
        x1 = tempSum * optF.x1Field / ((1 << optF.x1Size) - 1) + optF.getX1LowBound();
        //解码x2
        tempSum = 0;
        for (int j = optF.x1Size + optF.x2Size - 1, k = 0; j >= optF.x1Size; j--, k++)
        {
            tempSum += (1 << k) * individual[i][j];
        }
        x2 = tempSum * optF.x2Field / ((1 << optF.x2Size) - 1) + optF.getX2LowBound();
        //返回fit值
        return optF.evaluate(x1, x2);
    }

    //计算个体的适应力
    void evaluateFit(std::vector<std::vector<int>> &individual, double fit[])
    {
        for (int i = 0; i < N; i++)
        {
            //个体i的基因型解码计算适应力
            fit[i] = (float) 1.0 / decode(individual, i);
        }
    }

    //保留精英个体
    void saveElite(std::vector<std::vector<int>> &individual, int eliteNum[], std::vector<std::vector<int>> &eliteGene,
                   const double fit[], double eliteValue[])
    {
        //精英保留E个最优个体直接进入下一代种群
        for (int i = 0; i < E; i++)  //循环找出E个精英交叉池编号  i表示第i个精英
        {
            double maxFit = 0;
            for (int j = 0; j < N; j++)  //遍历所有个体适应值，找出第i个精英   j表示当前个体
            {
                if (maxFit < fit[j])   //判断当前最大适应值是否小于等于当前个体适应值
                {
                    bool repeat = false;
                    for (int k = 0; k < E; k++)  //遍历所有精英，判断当前精英是否重复
                    {
                        if (j == eliteNum[k])
                        {
                            repeat = true;
                            break;
                        }

                    }
                    if (!repeat)   //当前精英未重复
                    {
                        maxFit = fit[j];
                        eliteValue[i] = (float) 1.0 / maxFit;
                        eliteNum[i] = j;
                        for (int bit = 0; bit < optF.x1Size + optF.x2Size; bit++)
                        {
                            eliteGene[i][bit] = individual[j][bit];
                        }

                    } else continue;
                }
            }
        }
    }

    //把父代Individual[][]中第i个体放入池pool中，
    void copyIn(std::vector<std::vector<int>> &pool, std::vector<std::vector<int>> &individual, int i, int j)
    {
        for (int k = 0; k < optF.x1Size + optF.x2Size; k++)
            pool[i][k] = individual[j][k];
    }

    //建立交叉池
    void createPool(std::vector<std::vector<int>> &individual, std::vector<std::vector<int>> &pool, const double fit[])
    {
        //首先，计算个体的适应度占总适应度的比值pi，然后，计算个体的累计适应度qi。
        auto pi = new double[N];
        auto qi = new double[N];

        //求总适应度sum_fit
        double sumFit = std::accumulate(fit, fit + N, 0.0);
        //遍历所有个体适应度，求每个个体/比值pi，每个个体累计适应度qi
        for (int i = 0; i < N; i++)
        {
            pi[i] = fit[i] / sumFit;
            if (i == 0)
                qi[i] = pi[i];
            else
                qi[i] = qi[i - 1] + pi[i];
        }
        for (int i = 0; i < N; i++)
        {
            //轮盘赌选择概率p
            double p = u(e);
            for (int j = 0; j < N; j++)
                if (qi[j] > p)
                {
                    // 把N个体放入交叉池
                    copyIn(pool, individual, i, j);
                    break;
                }

        }
    }


    //两点交叉，单点变异
    void indiCrossVariation(std::vector<std::vector<int>> &pool)
    {
        for (int i = 0; i < N / 2; i+=2)
        {
            //随机值小于交叉概率则进行两点交叉
            if (this->Pc >= u(e))
            {
                //在交叉池内进行两点交叉
                unsigned left = e() % optF.x1Size + optF.x2Size, right = e() % optF.x1Size + optF.x2Size;
                if (left > right)
                    std::swap(left, right);
                for (auto j = left; j <= right; j++)
                    std::swap(pool[i][j], pool[i + 1][j]);
                //对交叉的个体的两个后代进行单点变异
                for (int j = 0; j < optF.x1Size + optF.x2Size; j++)
                {
                    if (this->Pm >= u(e))
                    {
                        pool[i][j] ^= 1;
                    }
                    if (this->Pm >= u(e))
                    {
                        pool[i+1][j] ^= 1;
                    }
                }
            }
        }
    }

    //计算交叉池个体的适应力
    void calculationFit(std::vector<std::vector<int>> &doubleIndi, double doubleFit[])
    {
        //计算适应力
        for (int i = 0; i < N * 2; i++)
        {
            doubleFit[i] = (float)1.0 / decode(doubleIndi, i);
        }
    }

    //计算累计适应度
    void calculationDoubleQi(double doubleFit[], double doublePi[], double doubleQi[])
    {
        //求总适应度sum_fit
        double sumFitAll = std::accumulate(doubleFit, doubleFit + 2 * N, 0.0);
        //遍历所有个体适应度，求每个个体/比值，每个个体累计适应度
        for (int i = 0; i < N * 2; i++)
        {
            doublePi[i] = doubleFit[i] / sumFitAll;
            if (i == 0)
                doubleQi[i] = doublePi[i];
            else
                doubleQi[i] = doubleQi[i - 1] + doublePi[i];
        }
    }

    void updatePopulation(std::vector<std::vector<int>> &individual, std::vector<std::vector<int>> &pool,
                          std::vector<std::vector<int>> &eliteGene)
    {
        //把父代和子代放入种群池doubleIndi
        std::vector<std::vector<int>> doubleIndi(2 * this->N, std::vector<int>(optF.x1Size + optF.x2Size));
        int i, j;
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < optF.x1Size + optF.x2Size; j++)
            {
                doubleIndi[i][j] = individual[i][j];
            }
        }
        for (i = N; i < 2 * N; i++)
            for (j = 0; j < optF.x1Size + optF.x2Size; j++)
            {
                doubleIndi[i][j] = pool[i - N][j];
            }
        //创建适应度数组
        auto doubleFit = new double[2 * N];
        //创建个体的适应度占总适应度的比值
        auto doublePi = new double[2 * N];
        //创建累计适应度数组
        auto doubleQi = new double[2 * N];
        //计算所有个体适应度,返回 doubleFit
        calculationFit(doubleIndi, doubleFit);
        //计算累计适应度，返回DoubleQi
        calculationDoubleQi(doubleFit, doublePi, doubleQi);
        //从种群 DoubleIndi中轮盘赌选择 N-E 个个体进入下一代 individual[N][x1_size+x2_size]
        for (i = 0; i < N - E; i++)
        {
            double p2 = u(e);
            for (j = 0; j < N * 2; j++)
            {
                if (doubleQi[j] > p2)
                {
                    copyIn(individual, doubleIndi, i, j);
                    break;
                }
            }
        }
        //精英个体EliteGene[E][x1_size+x2_size]直接加入下一代 individual[N][x1_size+x2_size]
        for (j = 0; j < E; j++, i++)
        {
            for (int k = 0; k < optF.x1Size + optF.x2Size; k++)
                individual[i][k] = eliteGene[j][k];
        }
        delete[]doubleFit;
        delete[]doublePi;
        delete[]doubleQi;
    }

    void decodeX(std::vector<std::vector<int>> &eliteGene, double *x, int l, int r, int xn)
    {
        for (int i = 0; i < E; i++)
        {
            int tempSum = 0;
            for (int j = r, k = 0; j >= l; j--, k++)
            {
                tempSum += (int) pow(2, k) * eliteGene[i][j];
            }
            if (xn == 1)
            {
                x[i] = tempSum * optF.x1Field / ((1 << (r - l + 1)) - 1) + optF.getX1LowBound();
            } else if (xn == 2)
            {
                x[i] = tempSum * optF.x2Field / ((1 << (r - l + 1)) - 1) + optF.getX2LowBound();
            }

        }

    }

//输出精英的适应值和基因码
    void outputElite(double eliteValue[], std::vector<std::vector<int>> &eliteGene, int current)
    {
        auto x1 = new double[N], x2 = new double[N];
        decodeX(eliteGene, x1, 0, optF.x1Size - 1, 1);
        decodeX(eliteGene, x2, optF.x1Size, optF.x1Size + optF.x2Size - 1, 2);
        //输出精英的适应值
        for (int i = 0; i < E; i++)
        {
            std::cout << "种群第 " << current << " 代第" << i + 1 << "个精英个体函数值：" << eliteValue[i];
            std::cout << " | (x1,x2)=(" << x1[i] << "," << x2[i] << ")" << " | 基因型：";
            for (int j = 0; j < optF.x1Size + optF.x2Size; j++)
            {
                std::cout << eliteGene[i][j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

public:
    OptFunction optF;

    GeneticAlgOptimizer(int n, int eliteN, int times, double Pc, double Pm, double ep)
    {
        this->N = n; //种群数量
        this->E = eliteN; // 精英数量
        this->Times = times; //迭代次数
        this->Pc = Pc; //两点交叉概率
        this->Pm = Pm; //单点变异概率
        this->epsilon = ep; //精度
    }

    void Optimize(std::function<double(double, double)> f, double x1Low, double x1High, double x2Low, double x2High)
    {
        //定义函数以及函数的上下界
        OptFunction optF(f, x1Low, x1High, x2Low, x2High, this->epsilon);
        this->optF = optF;
        //创建N个初始个体
        std::vector<std::vector<int>> individual(this->N, std::vector<int>(optF.x1Size + optF.x2Size));
        //存储每个个体的适应力
        auto fit = new double[this->N];
        //创建E个精英个体
        std::vector<std::vector<int>> eliteGene(this->E, std::vector<int>(optF.x1Size + optF.x2Size));
        //精英个体编号
        auto eliteNum = new int[this->E];
        //精英个体适应值
        auto eliteValue = new double[this->E]();
        //随机创建N个个体的初始基因型
        for (int i = 0; i < this->N; i++)  //遍历每个个体v
        {
            for (int j = 0; j < optF.x1Size + optF.x2Size; j++)    //每个个体随机x1_size+x2_size位01基因序列
            {
                individual[i][j] = (e() % 2);
            }
        }
        //迭代
        for (int current = 0; current < Times; current++)
        {
            //计算个体的适应力
            evaluateFit(individual, fit);
            //保留E个适应度最高的精英：double Elite[E] ，以及它的基因型
            saveElite(individual, eliteNum, eliteGene, fit, eliteValue);
            //轮盘赌选择创建交叉池
            std::vector<std::vector<int>> pool(N, std::vector<int>(optF.x1Size + optF.x2Size));
            createPool(individual, pool, fit);
            //种群个体间基因在交叉池两点交叉以及单点变异 产生在交叉池中的子代
            indiCrossVariation(pool);
            //把父代和子代合并，进行轮盘赌选择N-E个个体，更新种群个体
            updatePopulation(individual, pool, eliteGene);
            //输出当前代的精英表现形（适应值）以及基因型
            outputElite(eliteValue, eliteGene, current);
        }

    }

};


int main()
{
    GeneticAlgOptimizer geneticAlgOptimizer(1000, 3, 100, 0.9, 0.06, 0.0001);
    std::cout << "# 求解函数F1优化问题 #" << std::endl;
    geneticAlgOptimizer.Optimize(f1, -100, 100, -100, 100);
    std::cout << "# 求解函数F2优化问题 #" << std::endl;
    geneticAlgOptimizer.Optimize(f2, -10, 10, -10, 10);
    std::cout << "# 求解函数F3优化问题 #" << std::endl;
    geneticAlgOptimizer.Optimize(f3, -100, 100, -100, 100);
//    函数4不能使用函数值的倒数作为适应度，不考虑其优化问题
//    std::cout << "# 求解函数F4优化问题 #" << std::endl;
//    geneticAlgOptimizer.Optimize(f4, -500, 500, -500, 500);
    std::cout << "# 求解函数F5优化问题 #" << std::endl;
    geneticAlgOptimizer.Optimize(f5, -5.12, 5.12, -5.12, 5.12);
    std::cout << "# 求解函数F6优化问题 #" << std::endl;
    geneticAlgOptimizer.Optimize(f6, -32, 32, -32, 32);
    return 0;
}
