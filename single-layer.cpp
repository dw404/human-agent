#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include<string.h>
/* Random function */
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df          /* constant vector a */
#define UPPER_MASK 0x80000000        /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff        /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define NAMEOUT     "K4b075r5Q2"
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)
static unsigned long mt[NN];        /* the array for the state vector  */
static int mti = NN + 1;            /* mti==NN+1 means mt[NN] is not initialized */
int lc, ld, pc, pd;/*lc为选择合作的内在学习者
	                 pc为选择合作的外部模仿者
	 				 ld为选择背叛的内在学习者
					 pd为选择背叛的外部模仿者*/
void sgenrand(unsigned long seed) {
	int i;
	for (i = 0; i < NN; i++) {
		mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
		mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
	}
	mti = NN;
}
void lsgenrand(unsigned long seed_array[]) {
	int i; for (i = 0; i < NN; i++) mt[i] = seed_array[i]; mti = NN;
}
double genrand() {
	unsigned long y;
	static unsigned long mag01[2] = { 0x0, MATRIX_A };
	if (mti >= NN)
	{
		int kk;
		if (mti == NN + 1) sgenrand(4357);
		for (kk = 0; kk < NN - MM; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + MM] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (; kk < NN - 1; kk++) {
			y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
			mt[kk] = mt[kk + (MM - NN)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[NN - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[NN - 1] = mt[MM - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
	y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
	return y;
}
double randf() { return ((double)genrand() * 2.3283064370807974e-10); } //返回随机0~1的小数
long randi(unsigned long LIM) { return((unsigned long)genrand() % LIM); } //返回不大于LIM的随机整数
/* End of random */

#define RANDOMIZE 3145215
#define L 200 //L是边界大小
#define K 0.1 //K为噪声大小，此处就定为0.1
#define NEINUM 4 //有四个邻居
#define SIZE (L*L)
#define MC_STEPS 200000//一共迭代多少次
#define REC_STEPS 5000  //最后5000个的平均值作为实验结果
#define REFRESH_PRE 100	 /* Control the frequency of refresh screen, 1 for all the time */
int BM_NUM[SIZE]; //外部模仿者周围的BM数量,若它本身就是BM，则设定BM_NUM[x]为5
double b; //b是选择背叛的收益，也就是社会困境的强度
int net[SIZE][4];	/*net[x][0~3]代表x的四个邻居的位置*/ 
int cooperator, defector;  
char frequency[100];//字符数组用来给txt文件命名
double β = 1; 
double A = 0.5; //预期收入
double μ = 0.5;  //BM模型内在学习者的密度
double t = 0.2;
double p[SIZE];//第i个玩家（最多有size个）选择合作的概率
struct Strategy 
{
	int s;  //s为游戏策略，0合作，1背叛
	int type;  //type为更新策略，0为内在学习者BM，1为外部模仿者,0内1外
}stra[SIZE + 1];

/* Payoff matrix and its update */
double payoff_matrix[2][2] = { {1,0},
							   {1,0} };
/*收益矩阵 payoff_matrix[我的策略][你的策略]
0合作，1背叛
    C                D
C  [0][0]=1       [0][1]=0
  （都合作）    （我合作你背叛）
D  [1][0]=b       [1][1]=0
（我背叛你合作）   （都背叛）*/

/*初始先把收益矩阵的合作-背叛收益定为1，之后根据需要定为b*/
#define update_matrix(b) payoff_matrix[1][0]=b;  

//找到player x（第i+1行，j+1列）周围的四个邻居的位置
void prod_neighbors()
 {
	int i, j, x;
	for (i = 0; i < L; i++) {
		for (j = 0; j < L; j++) {
			x = i * L + j;//（从0开始）在网格内从左往右，从上往下的第几个
			net[x][0] = i * L + ((j - 1 + L) % L);	/*left*/
			net[x][1] = ((i - 1 + L) % L) * L + j;	/*up*/
			net[x][2] = ((i + 1) % L) * L + j;	/*down*/
			net[x][3] = i * L + ((j + 1) % L);	/*right*/
		}
	}
}


void picture(int step)//以点阵图的形式记录当前的各个玩家的游戏策略
{
    sprintf(frequency, "step=%d,distribution.txt",step);
	FILE* p3 = fopen(frequency, "w");
				for (int i = 0; i < L; i++)
				{
					for (int j = 0; j < L; j++)
					{
						if (stra[i * L + j].s == 0)//第i+1行第j+1列的玩家
						{
							if (stra[i * L + j].type == 0)
							{
								fprintf(p3, "0"); //选择合作的内在学习者标记为0
							}
							else
							{
								fprintf(p3, "1"); //选择合作的外部模仿者标记为1
							}
						}
						else
						{
							if (stra[i * L + j].type == 0)
							{
								fprintf(p3, "2"); //选择背叛的内在学习者标记为2
							}
							else
							{
								fprintf(p3, "3"); //选择背叛的外部模仿者标记为3
							}
						}
						if (j != L - 1)//如果不是最后一行，每输出一个数就空一格
						{
							fprintf(p3, " ");
						}
			        }
				fprintf(p3, "\n");
                }
	fclose(p3);
}

//对每个player进行策略以及位置的初始化
void init()
{
	cooperator = defector = 0;
	for (int i = 0; i < SIZE; i++)
		p[i] = 0.5;
	
	for (int i = 0; i < SIZE; i++)//完全随机的生成方式
	{
		//假如μ=0.3，那么randf()随机生成的0~1之间的小数小于0.3的概率就是0.3，那么内在学习者生成者所占比例就是0.3
		if (randf() < μ)
		{
			stra[i].type = 0;//内在学习者BM策略
		}
		else
		{
			stra[i].type = 1;//外部模仿者策略
		}
		//初始的合作率定为0.5
		if (randf() < 0.5)
		{
			stra[i].s = 0;//合作
		}
		else
		{
			stra[i].s = 1;//背叛
		}
		if (stra[i].s == 0)
			cooperator++;
		else
			defector++;
	}
    //picture(0);
}

void BM_num()//找到每个玩家周围有几个BM模型，若它自己就是BM模型，则为5；
{
	int j;
	for(j = 0; j < SIZE; j++)//遍历每个节点
	{
	    if(stra[j].type==0)//如果是BM模型的内在学习者
	        BM_NUM[j]=5;//人为设定为5，不代表实际意义
	    else//如果是外部模仿者
	   {
		    int i;
		    for(i = 0; i < NEINUM - 1; i++)//遍历四个邻居
		    {
			if(stra[net[j][i]].type==0)//若邻居是BM
			BM_NUM[j]++;//BM数量加1
		    }
	   }
	}
}

//与自己的邻居进行对比，算出收益
double payoff(int x)
{
	int i;
	double pay = 0;
	if(stra[x].type==0)//如果是自己是内在学习者BM模型，正常算
	{
		for (i = 0; i < NEINUM; i++)
	        pay += payoff_matrix[stra[x].s][stra[net[x][i]].s];
	//stra[x].s是自己的游戏策略，net[x][i]记录x的四个邻居的位置，stra[net[x][i]].s记录x四个邻居的游戏策略
	}
	else if(stra[x].type==1)//如果自己是外部模仿者
	{
		if(BM_NUM[x]==0)//如果周围四个人都是外部模仿者，还是正常算
		{
			for (i = 0; i < NEINUM; i++)
	            pay += payoff_matrix[stra[x].s][stra[net[x][i]].s];
		}
		else//如果周围有BM，那就只和内在学习者BM模型进行交互
		{
			for (i = 0; i < NEINUM; i++)
			{
				if(stra[net[x][i]].type==0)
				pay += payoff_matrix[stra[x].s][stra[net[x][i]].s];
			}
		}
	}
	return pay;
}


double stimu(int x)
{
	double s, r;//r为周围玩家的平均收益
	if(stra[x].type==0)
	r = payoff(x) / 4; //BM和周围四个玩家都交互，所以除以4
	else if(stra[x].type==1)
	{
		if(BM_NUM[x]==0) //普通玩家周围没有BM,只有普通玩家，还是正常算
		r = payoff(x) / 4;
		else
		r = payoff(x) / BM_NUM[x];
	}
	A = 0.5 * A + 0.5 * r;
	s = tanh(β * (r - A));
	return s;
}

void update_stra() 
{
	int i, j,player;
	double pay1, pay2;
	for (i = 0; i < SIZE; i++)
	{
		player = i;
		if (stra[player].type == 0) //内在学习者的更新策略
		{
			if (randf() <= p[player])//非常随机的改变，随机性的社会学习，所以用randf生成随机数
			{
				stra[player].s = 0;
			}
			else
			{
				stra[player].s = 1;
			}
		}
		else //外部模仿者的更新策略
		{
		    pay1 = payoff(player);  //自己的收益
			if(BM_NUM[player] == 0) //如果他周围没有BM，那就正常算
			{
				j = net[player][(int)randi(NEINUM)];  //随机生成player的一个邻居
			    pay2 = payoff(j);  //邻居的收益
			    if (stra[player].s != stra[j].s)
				{
					if (randf() < 1 / (1 + exp((pay1 - pay2) / K)))//也是非常随机的改变，随机性的社会学习，所以用randf生成随机数
					{
					    stra[player].s = stra[j].s;  //别人的好，学习别人的游戏策略
				    }
			    }
			}
			else //如果周围有BM，则只和BM进行交互
			{
				while (1)//无限循环，直到找到一个BM邻居后退出循环
				{
					j = net[player][(int)randi(NEINUM)];//随机找一个邻居
				    if(stra[j].type==0)//如果这个邻居是BM，则以他为对象进行模仿
				    {
						pay2 = payoff(j);
				        if (stra[player].s != stra[j].s)
				        {
							if (randf() < 1 / (1 + exp((pay1 - pay2) / K)))//也是非常随机的改变，随机性的社会学习，所以用randf生成随机数
					        {
						        stra[player].s = stra[j].s;  //别人的好，学习别人的游戏策略
					        }
				        }
				        break;
				    }
			    }

			}
			
		}
	
	}
	return;
}

/*
//普通交互模式
double payoff(int x)
{
	int i;
	double pay = 0;
	for (i = 0; i < NEINUM; i++)
		pay += payoff_matrix[stra[x].s][stra[net[x][i]].s];	
//stra[x].s是自己的游戏策略，net[x][i]记录x的四个邻居的位置，stra[net[x][i]].s记录x四个邻居的游戏策略
	return pay;
}

double stimu(int x)
{
	double s, r;
	r = payoff(x) / 4;
	s = tanh(β * (r - A));
	return s;
}


//更新游戏策略
void update_stra() 
{
	int i, j, player;
	double pay1, pay2;
	for (i = 0; i < SIZE; i++)
	{
		//player = (int)randi(SIZE);//此处有小问题，本来是想游戏顺序也随机，但是并没有保证每个人都能更新到
		player = i;
		if (stra[player].type == 0) //内在学习者的更新策略
		{
			if (randf() <= p[player])//非常随机的改变，随机性的社会学习，所以用randf生成随机数
			{

				stra[player].s = 0;
			}
			else
			{
				stra[player].s = 1;
			}
		}
		else //外部模仿者的更新策略
		{
			pay1 = payoff(player);  //自己的收益
			j = net[player][(int)randi(NEINUM)];  //随机生成player的一个邻居
			pay2 = payoff(j);  //邻居的收益
			if (stra[player].s != stra[j].s)
			{
				if (randf() < 1 / (1 + exp((pay1 - pay2) / K)))//也是非常随机的改变，随机性的社会学习，所以用randf生成随机数
				{
					stra[player].s = stra[j].s;  //别人的好，学习别人的游戏策略
				}
			}
		}
	
	}
	return;
}
*/

/*内在学习者BM模型的选择合作的概率*/
void calcul()
{
	int i, j;
	double s;
	for (i = 0; i < SIZE; i++)
	{
		s = stimu(i);
		if ((s >= 0) && stra[i].s == 0)
		{
			p[i] = p[i] + (1 - p[i]) * s;
		}
		else if ((s < 0) && stra[i].s == 0)
		{
			p[i] = p[i] + p[i] * s;
		}
		else if ((s >= 0) && stra[i].s == 1)
		{
			p[i] = p[i] - p[i] * s;
		}
		else
		{
			p[i] = p[i] - (1 - p[i]) * s;
		}
	}
}

//更新四种玩家的数量，并计算合作者与背叛者各自的总数量
void update_data() 
{
	int i;
	cooperator = defector = 0;
	lc = pc = ld = pd = 0;
	for (i = 0; i < SIZE; i++) {
		if (stra[i].s == 0)
		{
			cooperator++;//总的合作者
			if (stra[i].type == 0)
			{
				lc++;//选择合作的内在学习者
			}
			else
			{
				pc++;//选择合作的外部模仿者
			}
		}
		else
		{
			defector++;//总的背叛者
			if (stra[i].type == 0)
			{
				ld++;//选择背叛的内在学习者
			}
			else
			{
				pd++;//选择背叛的外部模仿者
			}
		}
	}
}

int main()
{

	sgenrand(RANDOMIZE);
	prod_neighbors();
	int x, step;
	double fc, afc = 0;  //fc为合作者的密度也就是合作率，afc为平均合作率
	double ifc, efc = 0;  //ifc为BM模型内在学习者的合作率,efc为外部模仿者的合作率
    b=1.025;
    //μ=0;
	/*μ是全局变量，通过μ的递增，可以探索不同的内在学习者的密度对于合作率的影响*/
	//for(b = 1.1; b <= 2.01; b = b + 0.1)
	{
		for (μ = 0; μ <= 1.01; μ = μ + 0.1)
			{
				init();//根据μ的不同，初始化不同的初始分布
				update_matrix(b);
				for (step = 1; step < MC_STEPS; step++)//根据之前定好的step大小开始循环
				{
					calcul();
					update_stra();
					update_data();
					fc = (double)cooperator / SIZE; //合作者的密度，总合作率
					/*if(μ==0)//没有BM，全是外部模仿者
					{
						ifc=0;
						efc = (double)pc / SIZE; //外部模仿者的合作率
					}
					else if(μ==1)//没有外部模仿者，全是BM
					{
						efc=0;
						ifc = (double)lc / SIZE; //BM的合作率
					}
					else//两者都有
					{
                        efc = (double)pc / (SIZE*(1-μ)); //外部模仿者的合作率
						ifc = (double)lc / (SIZE*μ); //BM的合作率
					}
					*/
					if (step > MC_STEPS - REC_STEPS - 1) //排除干扰,只考虑最后五千步
						afc += fc; //计算最后五千步的合作率之和
					/*if (step <=100) 
					{
						sprintf(frequency, "b=2, μ=%lf, step new.txt",μ);
						FILE* p7 = fopen(frequency, "a");
						fprintf(p7, " %lf%% \n",  fc*100);
						fclose(p7);
					}*/
					if (step % REFRESH_PRE == 0)  //如果是步数是一百的倍数，就打印步数与合作率
					{
						printf("\rStep=%d  \tC=%lf%%  μ=%lf	", step, fc * 100, μ);
				        /*sprintf(frequency, "b=2, μ=%lf, step new.txt",μ);
				        FILE* p7 = fopen(frequency, "a");
			            fprintf(p7, " %lf%%  \n",  fc*100);
						fclose(p7);*/
					}
					if ((cooperator==0) || (defector==0)) //当只剩下C或只剩下D时
					{		
						if (step++ < MC_STEPS - REC_STEPS)//若程序还未进行到最后一步
						/*如果只剩背叛者，那么合作者数量为0，合作率afc为0
						如果只剩合作者，那么合作者数量不为零，合作率为1，afc也为1*/
							afc = cooperator ? 1 : 0;
							/*我对于上一句的修改：
							{
								if(cooperator==0)
								afc=0;
								else if(defector==0)
								afc=1;
							}*/
						break;//退出循环
					}
					/*if(step == 1)
					picture(1);
					if(step == 10)
					picture(10);
					if(step == 30)
					picture(30);
					if(step == 100)
					picture(100);
					if(step ==  MC_STEPS - 2)
					picture(MC_STEPS - 2);*/
				}//step循环结束
				afc /= step + REC_STEPS - MC_STEPS;//计算最后五千步的平均合作率
				printf("\rb=%lf \tavg_C=%lf%% μ=%g", b, afc * 100, μ);
				sprintf(frequency, "dynamic A, b=1.025, μ=0-1, afc new.txt");
				FILE* p7 = fopen(frequency, "a");
			    fprintf(p7, " %lf %lf %lf%%  \n",  b, μ, afc*100);
				//fprintf(p7, " %lf %lf %lf \n",  μ, ifc, efc);
				fclose(p7);
				printf("\n");
               //_fcloseall();
			}//μ循环结束

    }//b循环结束
	
	return 0;
}
