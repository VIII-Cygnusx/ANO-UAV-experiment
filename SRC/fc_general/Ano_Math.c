/******************** +++++++++++++++++++++++++++ ********************************
 * 描述    ：快速计算
**********************************************************************************/

#include "Ano_Math.h"


float my_abs(float f)
{
	if (f >= 0.0f)
	{
		return f;
	}

	return -f;
}

REAL fast_atan2(REAL y, REAL x) 
{
	#define ATANmagik1 -0.05030176426f                                          // Double: -0.05030176425872175099
    #define ATANmagik2 -6.988836621f                                            // Double: -6.9888366207752135
    #define ATANmagik3  3.145599955e-7f                                         // Double:  3.14559995508649281e-7
    #define ATANmagik4  2.844463688f                                            // Double:  2.84446368839622429
    #define ATANmagik5  0.8263997833f                                           // Double:  0.826399783297673451
    #define ATANmagik6  0.1471039134f                                           // Double:  0.1471039133652469065841349249
    #define ATANmagik7  0.6444640677f                                           // Double:  0.644464067689154755092299698
    float res, rsq, absX, absY;
    absX = ABS(x);
    absY = ABS(y);
    res  = _MAX(absX, absY);
    if (res) res = _MIN(absX, absY) / res;
    else res = 0.0f;
    rsq = res * res;
    res = ATANmagik1 * (ATANmagik2 + res) * (ATANmagik3 + res) * (ATANmagik4 + ATANmagik5 * res + rsq) /
          (1.0f + ATANmagik6 * res + ATANmagik7 * rsq);
    if (absY > absX) res = MY_PPPIII_HALF - res;
    if (x < 0) res = MY_PPPIII - res;
    if (y < 0) res = -res;
    return res;
}

float my_atan(float x, float y)
{
	return fast_atan2(y, x);
}

////计算浮点数平方
//float my_pow(float a)
//{
//	return a*a;
//}

float my_sqrt_reciprocal(float number)
{

	long i;
	float x, y;

	x = number * 0.5F;
	y = number;
	i = * ( long * ) &y;
	i = 0x5f3759df - ( i >> 1 );

	y = * ( float * ) &i;
	y = y * ( 1.5f - ( x * y * y ) );
	y = y * ( 1.5f - ( x * y * y ) );
	
	return y;
}

//快速平方根算法
float my_sqrt(float number)
{

	return number * my_sqrt_reciprocal(number);
}

#define ONE_PI   (3.14159265)
#define TWO_PI   (2.0 * 3.14159265)
#define ANGLE_UNIT (TWO_PI/10.0)

double my_sin(double x)
{
	#define sinPolyCoef3 -1.666568107e-1f
	#define sinPolyCoef5  8.312366210e-3f
	#define sinPolyCoef7 -1.849218155e-4f
	#define sinPolyCoef9  2.600054768e-6f
    int32_t xint = x;
    if (xint < -32 || xint > 32) return 0.0f;                               // Stop here on error input (5 * 360 Deg)
    while (x >  MY_PPPIII) x -= (2.0f * MY_PPPIII);                                 // always wrap input angle to -PI..PI
    while (x < -MY_PPPIII) x += (2.0f * MY_PPPIII);
    if (x >  (MY_PPPIII_HALF)) x =  (MY_PPPIII_HALF) - (x - (MY_PPPIII_HALF));   // We just pick -90..+90 Degree
    else if (x < -(MY_PPPIII_HALF)) x = -(MY_PPPIII_HALF) - ((MY_PPPIII_HALF) + x);
    float x2 = x * x;
    return x + x * x2 * (sinPolyCoef3 + x2 * (sinPolyCoef5 + x2 * (sinPolyCoef7 + x2 * sinPolyCoef9)));
}

float my_cos(double rad)
{
	return my_sin(rad+MY_PPPIII_HALF);
}

float my_deadzone(float x,float ref,float zoom)
{
	float t;
	if(x>ref)
	{
		t = x - zoom;
		if(t<ref)
		{
			t = ref;
		}
	}
	else
	{
		t = x + zoom;
		if(t>ref)
		{
			t = ref;
		}
	}
  return (t);
}

float my_deadzone_2(float x,float ref,float zoom)
{
	float t;
	
	if( x> (-zoom + ref) && x < (zoom + ref) )
	{
		t = ref;
	}
	else
	{
		t = x;
	}

  return (t);
}

float my_HPF(float T,float hz,float x,float zoom,float range,float *zoom_adj)
{
// 	if( ABS(x-ref) <zoom )
// 	{
// 	 // *zoom_adj += hz *T *3.14f *(x - *zoom_adj);
// 		*zoom_adj += ( 1 / ( 1 + 1 / ( hz *6.28f *T ) ) ) *LIMIT( (x - *zoom_adj),(- zoom),(zoom) );
// 		*zoom_adj = LIMIT(*zoom_adj,-zoom,zoom);
// 	}
	if( ABS(x ) < 0.5f *range *zoom )
	{
		hz *= 1.2f; 
	}
	else if( ABS(x ) < range *zoom )
 	{
 		hz *= 0.8f; 
 		
	}
	else if( ABS(x ) < zoom )
	{
		hz *= 0.5f; 
	}
	else if( ABS(x ) < 2*zoom )
	{
		hz *= 0.2f; 
	}
	else
	{
		hz *= 0.1f; 
	}
	
	*zoom_adj += ( 1 / ( 1 + 1 / ( hz *6.28f *T ) ) ) *(x - *zoom_adj);
	*zoom_adj = LIMIT(*zoom_adj,-range *zoom,range *zoom);
	return (x - *zoom_adj);

}


double To_180_degrees_db(double x)
{
	return (x>180?(x-360):(x<-180?(x+360):x));
}

void length_limit(float *in1,float *in2,float limit,float out[2])
{
	float l = my_2_norm(*in1,*in2);
	float l_lim = LIMIT(l,0,limit);
	
	if(l==0)
	{
		out[0] = out[1] = 0;
	}
	else
	{
		out[0] = l_lim/l *(*in1);
		out[1] = l_lim/l *(*in2);
	}
}

float fifo(u8 arr_num,u8 *cnt,float *arr,float in)
{
	*(arr + *cnt) = in;
	*cnt += 1;
	if(*cnt>=arr_num) *cnt = 0;
	
	return (*(arr + *cnt));
}

//=========================================
//====================vector===============
/*
|x2|    |cosx,-sinx|   |x1|
|  | =  |          |   |  |
|y2|    |sinx, cosx|   |y2|
*/
void rot_vec_2(float in[2],float sinx,float out[2]) //x = +-90度旋转，取sin x
{
	out[0] = in[0] *my_sqrt(1-my_pow(sinx)) - in[1] *sinx;
	out[1] = in[1] *my_sqrt(1-my_pow(sinx)) + in[0] *sinx;
}

/*
va x vb = |va||vb| *sinx *vn,vn为垂直于va,vb的单位向量
归一化计算，取 va x vb = sinx;
*/

float vec_2_cross_product(float in1[2],float in2[2]) //正负为in1->in2 夹角逆时针
{
	return (in1[0] *in2[1] - in1[1] *in2[0]);
}

float vec_2_dot_product(float in1[2],float in2[2]) //正负为in1->in2 夹角（空间实际夹角）
{
	return (in1[0]*in2[0] + in1[1]*in2[1]);
}

/*
A x B = (AyBz - AzBy)i + (AzBx - AxBz)j + (AxBy - AyBx)k
			
*/
	
void vec_3_cross_product_err_sinx(float in1[3],float in2[3],float out[3]) //输出xyz误差夹角x 的sin(x)，右手螺旋
{
	out[0] = (in1[1] * in2[2] - in1[2] * in2[1] ) ;
	out[1] = (in1[2] * in2[0] - in1[0] * in2[2] ) ;
	out[2] = (in1[0] * in2[1] - in1[1] * in2[0] ) ;
}

float vec_3_dot_product(float in1[3],float in2[3]) //正负为in1->in2 夹角（空间实际夹角）
{
	return (in1[0]*in2[0] + in1[1]*in2[1] + in1[2]*in2[2]);
}

//向量乘矩阵转置
void Vec3f_Mul_MatrixT(float vec_in[3],float mat_in[3][3],float vec_out[3])
{
	//
	for(u8 i = 0;i<3;i++)
	{
		float temp = 0;
		for(u8 j = 0;j<3;j++)
		{
			
			temp += vec_in[j] *mat_in[i][j];
		}
		vec_out[i] = temp;
	}
}

