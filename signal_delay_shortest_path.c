#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#define ARRAYSIZE(x)  (sizeof(x)/sizeof(*(x)))
#define N 10000;
#define V 1658
float afunction(int a);
float bfunction(int b);
float max(float a, float b);
int sfunction(int s);
float min(float a, float b);
void dijkstra(float graph[V][V], int src, float Tx[5000], int sink);


float graph[V][V];
int main(void)
{
   /* Initially using an excel macro, all cell values of excel are multiplied by N and TRUNC function is used
   This is done to avoid the decimal points to be confused as delimiters */
   int max_row;
   const char filename[] = "cheninput.csv";
   /* Open the data file. */
   FILE *file = fopen(filename, "r");
   float array_read[5000][23];
   if ( file )
   {
      int array[5000][22];
      
      size_t i, j, k;
      char buffer[BUFSIZ], *ptr;
      /* Read each line from the file using fgets */
      for ( i = 0; fgets(buffer, sizeof buffer, file); ++i )
      {
         /*
           Parse the comma-separated values from each line into 'array'. */
         for ( j = 0, ptr = buffer; j < ARRAYSIZE(*array); ++j, ++ptr )
         {
            array[i][j] = (int)strtol(ptr, &ptr, 10);
         }
      }
      max_row=i;
      fclose(file);
     /* The array elements are fed into a new array which is float and contains the original values */
     
     
      for ( j = 0; j < i; ++j )
      {
         for ( k = 0; k < ARRAYSIZE(*array); ++k )
         {
            array_read[j][k]=(float)array[j][k]/N;
         }
      }      
   }
   else /* If the file does not exist, error message pops up. */
   {
      perror(filename);
   }

/* indegree and outdegree calculation begins: for unsignalised nodes */
int i,j;   
int nodedata[1658][3]={0};
for (i=0;i<=1657;i++)
{
	nodedata[i][0]=i+1;
	for(j=0;j<=max_row;j++)
	{
		if(array_read[j][1]==i)
		{
			nodedata[i][1]++; /*outdegree */
		}
        if(array_read[j][2]==i)
        {
			nodedata[i][2]++; /*indegree */
        }
    }
}


	/* code begins*/
                           
	/* shortest path algorithm begins */
	int sourcenode=0;
  int sinknode=0;
	printf("\nEnter source node = ");
	scanf("%d",&sourcenode);
  printf("\nEnter sink node = ");
  scanf("%d",&sinknode);
  sinknode=sinknode-1;
  int realsrc;
  realsrc = sourcenode-1;
  int realsink;
  realsink=sinknode-1;
  //int realsrc = sourcenode-1;
/* computong and saving travel times in the array Tx, average speed values in the array Sx,
and penalties in the array Px */
float To, x, cap, rat1, alpha, beta, term1,term2,term0,Tt,Tef;
float Tx[5000] = {0.0};
float Sx[5000] = {0.0};
float Te[5000] = {0.0};
float So,St;
float rsq2c, rat2;
int f,s,num_lane;
float a,b,d1x,d1,d2x,d2;
float one = 1;
float ptone = 0.01;
int k;
for(k=0;k<4224;k++)
  {
    
    To=array_read[k][6];
    So=array_read[k][8];
    x=array_read[k][21];
    cap=array_read[k][7];
    rsq2c=array_read[k][14];
    num_lane=(int)array_read[k][11];
    s=sfunction(num_lane);
    a=afunction(num_lane);
    b=bfunction(num_lane);
    rat1=(float) (x/cap);
    rat2=(float) (x/s);
    term2= (1-rat2);
    term0= max(term2,ptone);
    d1=(rsq2c/term0);
    d1x=((a*d1)+b);
    d2=array_read[k][19];
    d2x=d2*min(rat2,one);
    alpha=array_read[k][9];
    beta=array_read[k][10];
    term1=((alpha*pow(rat1,beta))+1);
    Tx[k]=To*term1;
    Tt=Tx[k];
    Te[k]=(Tt+d1x+d2x);
    Tef=Te[k];
    Sx[k]=(So/term1);
    St=Sx[k];
    f=k+1;
    //printf("\ncost %d = %f\n",f,Tef);
    //printf("travel time %d = %f\n",f,Tt);
    //printf("%d d1x\t %f \td2x\t %f\t\n",f,d1x,d2x);
    //printf("average speed %d = %f\n",f,St);

  }
int g,h;
  for(g= 0; g < V; g++)
  {
    for(h = 0 ; h < V; h++)
    {
      graph[g][h] = INT_MAX;
    }
  }
  
  for(i = 0; i < 4224; i++)
  {
    
    graph[(int)array_read[i][1]-1][(int)array_read[i][2]-1] = Te[i]; // storing the effective travel times between all pair of nodes
    graph[(int)array_read[i][2]-1][(int)array_read[i][1]-1] = Te[i];

    //printf("%f\n",graph[6][7]);
  }

  

	dijkstra(graph, realsrc, Te, sinknode);
 
	
	return 0 ;
}

float minDistance(float dist[], int sptSet[])
{
   int min_index, v;
   float min = INT_MAX;
 
   for (v = 0; v < V; v++)
     if (sptSet[v] == 0 && dist[v] <= min)
         min = dist[v], min_index = v;
 
   return min_index;
}

void dijkstra(float graph[V][V], int src, float Te[5000], int sink)
{
     float dist[V];
     int sptSet[V]; 
 	 int count, i, v, j;
   int pred[V];
 	 	
 
	   // Initialize all distances as INFINITE and stpSet[] to zero
     for ( i = 0; i < V; i++)
     {
        dist[i] = INT_MAX, sptSet[i] = 0;
     }
        // Distance of source vertex from itself is always 0
     dist[src] = 0;
     // Find shortest path for all vertices
     for (count = 0; count < V; count++)
     {
      /* Pick the minimum distance vertex from the set of vertices not
        yet processed. u is always equal to src in first iteration.*/
       int u = minDistance(dist, sptSet);
       // Mark the picked vertex as processed
       sptSet[u] = 1;
        
         // Update dist value of the adjacent vertices of the picked vertex.
       for (v = 0; v < V; v++)
        /* Update dist[v] only if is not in sptSet, there is an edge from 
         u to v, and total weight of path from src to  v through u is 
         smaller than current value of dist[v]*/
       {
        if (sptSet[v] != 1 && graph[u][v] != INT_MAX && dist[u] != INT_MAX  && dist[u]+graph[u][v] < dist[v])
        {
          dist[v] = dist[u] + graph[u][v];
          pred[v] = u; // storing the predecessor
        }
       }
     }
     int q;
   // printf("Vertex   Distance from Source  \t\tpred\n");
   for (i = 0; i < V; i++)
   {
      j=i+1;
      q=pred[i]+1;
      //printf("%d \t\t %f \t\t %d\n", j, dist[i],q);

   }
   

   //printing the path
  printf("\nThe shortest path is\n");
  int p;
  int pro[V];
  i=1;
  p = sink;
  while(p!=src)
  {    
    pro[i]=p+1;
    p=pred[p];
    i++;
  }
  int dest = pro[i]-1;
  int relsrc=src+1;
  if(relsrc!=0)
  printf("%d\n",relsrc);
  for(i=p;i>=1;i--)
  {
    if(i==1)
    {
      int dest = pro[i];
      if(dest!=0)
      printf("%d\n", dest);
    }
    else
    {
      if(pro[i]!=0)
      printf("%d\n",pro[i]);  
    }
  }
  system("pause");
}
 
float afunction(int a)
  { 
    if (a==2)
    {
    return 1.079;
    }
    else if(a==3) return 0.977;
    else if(a==4) return 0.870;
    else if(a==6) return 1.336;
  }

int sfunction(int s)
  {
  if (s==2)
  {
    return 4000;
  }
  else if(s==3) return 6000;
  else if(s==4) return 8000;
  else if(s==6) return 12000;
  }

float bfunction(int b)
  {
      if (b==2)
    {
    return 20.99;
    }
    else if(b==3) return 20.37;
    else if(b==4) return 33.80;
    else if(b==6) return 16.80;
  }   

float max(float a, float b)
{
      if (a>b) return a;
      else return b;
}

float min(float a, float b)
{
  if(a>b) return b;
  else return a;
}

