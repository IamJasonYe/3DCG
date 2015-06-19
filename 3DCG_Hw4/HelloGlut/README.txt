

【新增结构体】：
ASCDModel结构体：用于保存读入的模型数据

【新增函数】：
以下函数省略具体参数：

建立观察者视野：生成EM和PM
void observer();

建立视窗：生成WVM及PM的AR
void viewport();

展示图像：
void display();

求出两向量点积
float dot_product();

求出两向量叉乘
vec4 cross_product();

使一个向量归一化
vec4 normalize();

读取一个模型的数据
void readModel(string filename);

用model_matrix在world space生成物体
void create_object();

最后将原来2d的矩阵和向量升级成了3d的矩阵和向量

