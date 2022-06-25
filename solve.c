/*                                                                */
/*    訪問順制約付きTSP用プログラム                                    */
/*    C code written by K. Ando and K. Sekitani (Shizuoka Univ.)  */
/*                                                                */
/*                                                                */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "list.h"
#define MAX_N 10000   // 点の数の最大値
#define INF 100000000 // 無限大の定義
#define EPSILON 0.00000001 //ε 小さい正の値
#define SWAP(a,b){int temp; temp=(a); (a)=(b); (b)=temp; }   

struct point {
  int x;
  int y;
};

double dist(struct point p, struct point q) { // pとq の間の距離を計算 
  return sqrt((p.x-q.x)*(p.x-q.x)+(p.y-q.y)*(p.y-q.y));
}

double tour_length(struct point p[MAX_N], int n, int tour[MAX_N]) {
  int i;
  double sum=0.0;
  for(i=0;i<n;i++) sum+=dist(p[tour[i]],p[tour[(i+1)%n]]);
  return sum;// 総距離が関数の戻り値
}

void read_tsp_data(char *filename, struct point p[MAX_N],int *np, int prec[MAX_N], int *mp) {
  FILE *fp;
  char buff[500];
  int i;

  if ((fp=fopen(filename,"rt")) == NULL) {// 指定ファイルを読み込み用に開く
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }   

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // PRECEDENCE_CONSTRAINTS:で始まる行に出会う
	&&(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)!=0)) ; // まで読み飛ばす. 
  if(strncmp("PRECEDENCE_CONSTRAINTS:",buff,23)==0)  {
    sscanf(buff+24,"%d",mp);
    for(i=0;i<*mp;i++) fscanf(fp,"%d ",&prec[i]);
  } else {
    fprintf(stderr,"Error: There is no precedence constraint in file %s.\n",filename);
    exit(EXIT_FAILURE);
  }

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // DIMENSION で始まる行に出会う
	&&(strncmp("DIMENSION",buff,9)!=0)) ; // まで読み飛ばす. 
  sscanf(buff,"DIMENSION: %d",np);           // 点の数 *np を読み込む

  while((fgets(buff,sizeof(buff),fp)!=NULL)   // NODE_COORD_SECTIONで始まる
	&&(strncmp("NODE_COORD_SECTION",buff,18)!=0)) ; // 行に出会うまで, 
                                                        // 読み飛ばす. 
  for(i=0;i<*np;i++) {                       // i=0 から i=(*np)-1まで
    if(fgets(buff,sizeof(buff),fp)!=NULL) 
      sscanf(buff,"%*d %d %d",&(p[i].x),&(p[i].y)); // i番目の点の座標を読み込む
  }                                 

  fclose(fp);
}

void nn(struct point p[MAX_N],int n,int tour[MAX_N],int m, int prec[MAX_N]){
  int i,j,r,nearest=-1;
  int a;  
  int visited[MAX_N]; // 都市iが訪問済みならば1そうでなければ0
  double min;
  int l, k=2;
  int er=0;
  for(r=0;r<n;r++) visited[r]=0; // 最初は全都市は未訪問
  tour[0]=prec[0];         // 最初に訪問する都市は 0 
  visited[tour[0]]=1;      // 都市0は訪問済み

  for(i=0;i<n-1;i++) {
    a = tour[i];
    //最後に訪問した都市 a == tour[i]から最短距離にある未訪問都市nearestを
    //見つける
    min = INF;  
    for(r=0;r<n;r++) { 
      er = 0;   
      for(l = k; l< m; l++) {
          if(prec[l] == r) er=1;
        }
      if(er == 1) continue;
      //if(prec[k-1] == r) k++;
      if(!visited[r] && dist(p[a],p[r])<min){
	      nearest=r; //都市tour[i]から暫定的に最短距離にある未訪問都市をnearestとする
	      min = dist(p[a],p[r]); // その暫定最短距離をminとする
      }
      
       
    }
    if(prec[k-1] == nearest) k++;
    tour[i+1]=nearest; // i+1 番目に訪問する都市を nearest にして, 
    visited[nearest]=1;// nearest を訪問済みとする. 
    //    printf("tour   :"); show_array(tour,n);
    //    printf("visited:"); show_array(visited,n);
  }
  
}

void ci(struct point p[MAX_N],int n,struct list* tour, int m, struct list* unvisited) {
  int a,b,r;
  struct cell* i;
  struct cell* j;
  struct cell* rr;
  struct cell* nearest;
  double min_dist;
  

  // a= 0 に最も近い点を探す
  while(unvisited->head->next != unvisited->tail) {
    min_dist=INF;
    for(i=unvisited->head->next;i->next!=NULL;i=i->next) {
      r=i->data;
      for(j=tour->head->next;j->next!=NULL;j=j->next) {
        a=j->data;
        b=j->next->data;
        if(j->next == tour->tail) b=tour->head->next->data;
        if (dist(p[a],p[r])+dist(p[b],p[r])-dist(p[a],p[b])<min_dist) {
          rr=i;
          nearest=j;
          min_dist=dist(p[a],p[r])+dist(p[b],p[r])-dist(p[a],p[b]);
        }
        
      }
    }
    insertAfter(nearest,rr->data);
    erase(rr);
  }

}

void ni(struct point p[MAX_N],int n,struct list* tour, int m, struct list* unvisited) {
  int a,b,c,r,d=0;
  struct cell* i;
  struct cell* j;
  struct cell* rr;
  struct cell* aa;
  struct cell* nearest; 
  double min_dist;

  aa=tour->head->next;
  // a= 0 に最も近い点を探す
  while(unvisited->head->next != unvisited->tail) {
    min_dist=INF;
    a=aa->data;
    if(aa->prev==tour->head) {
      b=tour->tail->prev->data;
    } else {
      b=aa->prev->data;
    }
    
    if(aa->next==tour->tail) {
      c=tour->head->next->data;
    } else {
      c=aa->next->data;
    }
    
    for(i=unvisited->head->next;i->next!=NULL;i=i->next) {
      r=i->data;
      if (dist(p[a],p[r])<min_dist) {
        rr=i;
        min_dist=dist(p[a],p[r]);
      }
    }
    if (dist(p[a],p[r])+dist(p[b],p[r])-dist(p[a],p[b])<dist(p[a],p[r])+dist(p[c],p[r])-dist(p[a],p[c])) {
      insertBefore(aa,rr->data);
      aa=aa->prev;
      if(aa->prev==tour->head) aa=tour->tail->prev;
    } else {
      insertAfter(aa,rr->data);
      aa=aa->next;
      if(aa->next==tour->tail) aa=tour->head->next;
    }
    erase(rr);
  }
}

void TwoOpt(struct point p[MAX_N], int n, int tour[MAX_N], int m, int prec[MAX_N]){
  int a,b,c,d;
  int i,j,k,l,g,h;
  
  int s_tour,n_tour,n_prec;
  for(i=0;i<n;i++) {
    if(prec[0] == tour[i]) {
      s_tour = i;
      n_prec = i+1;
      break;
    }
  }

  double max_dif;
  int re = 1, con = 1;
  while(con) {
    con = 0;
    for(i=1;i<=2;i++) {
      if(tour[s_tour+i]==prec[n_prec]) {
        s_tour = (s_tour+i)%n;
        n_prec = (n_prec + 1)%m;
        con = 1;
      }
    }
    //printf("%d \n", prec[n_prec]);
    if(0 == n_prec) break;
    if(con) continue;
    re=1;
      while(re) {
      max_dif = -1;
      re = 0;

      for(i=s_tour;tour[(i+2)%n] != prec[n_prec];i++) {
          j=i+1;
          n_tour = (i+3)%n;
          for(k=i+2;tour[k] != prec[n_prec];k++){
              l = (k+1)%n;
              a = tour[i]; b = tour[j];
              c = tour[k]; d = tour[l];
              if(dist(p[a],p[b])+ dist(p[c],p[d]) > dist(p[a],p[c])+ dist(p[b],p[d])) {               
                  if(max_dif < (dist(p[a],p[b])+ dist(p[c],p[d])) - (dist(p[a],p[c])+ dist(p[b],p[d]))) {   
                      max_dif = (dist(p[a],p[b])+ dist(p[c],p[d])) - (dist(p[a],p[c])+ dist(p[b],p[d]));
                      g = j;
                      h = k;    
                  }
                  re = 1;
              }
          
          }
      }
      if(re) {
          while(g<h) {
              SWAP(tour[g],tour[h]); 
              g++;
              h--; 
          }
      }
    }

    n_prec = (n_prec+1)%m;
    if(0!=n_prec) con=1;
    s_tour = n_tour;
  }
}

void write_tour_data(char *filename, int n, int tour[MAX_N]){
  FILE *fp; 
  int i;
 
 // 構築した巡回路をfilenameという名前のファイルに書き出すためにopen
  if((fp=fopen(filename,"wt"))==NULL){ 
    fprintf(stderr,"Error: File %s open failed.\n",filename);
    exit(EXIT_FAILURE);
  }
  fprintf(fp,"%d\n",n);
  for(i=0;i<n; i++){
   fprintf(fp,"%d ",tour[i]);
  }
  fprintf(fp,"\n");
  fclose(fp);
}

int main(int argc, char *argv[]) {
  int n;                   // 点の数 
  int m;                   // 順序制約に現れる点の数
  struct point p[MAX_N];   // 各点の座標を表す配列 
  int tour[MAX_N];   // 巡回路を表現する配列
  int prec[MAX_N];   // 順序制約を表現する配列
  int i;
  int j, mp;
  struct list tour2;
  struct list unvisited;

  initialize(&tour2);
  initialize(&unvisited);

// 点の数と各点の座標を1番目のコマンドライン引数で指定されたファイルから読み込む
  read_tsp_data(argv[1],p,&n,prec,&m);

  for(i=0;i<m;i++) {
    insertBefore(tour2.tail,prec[i]);
  }  

  printNumbers(&tour2);

  for(int i=0;i<n;i++){
    int con=1;
    for(j=0;j<m;j++) {
      if(i==prec[j]) con=0;      
    }
    if(con) insertBefore(unvisited.tail,i);
  }

  if(argc != 2) {
    fprintf(stderr,"Usage: %s <tsp_filename>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  // 点の数と各点の座標を1番目のコマンドライン引数で指定されたファイルから読み込む

  //順序制約の確認
  //for(i=0;i<m;i++) printf("%d\n",prec[i]);

  // 最近近傍法による巡回路構築
  //nn(p,n,tour,m,prec);
  //Ceapest Insertionによる巡回路構築
  //ci(p,n,&tour2,m,&unvisited);
  //Nearest Insertionによる巡回路構築
  ni(p,n,&tour2,m,&unvisited);
  j=0;
  for(struct cell* i=tour2.head->next;i!=NULL;i=i->next) {
        tour[j] = i->data;
        j++;
    }
  //printNumbers(&tour2);
  printf("%5.1lf\n",tour_length(p,n,tour)); 

  TwoOpt(p,n,tour,m,prec);

  // ファイルに出力
  write_tour_data("tour1.dat",n,tour);
  printf("%5.1lf\n",tour_length(p,n,tour)); // 巡回路tourの長さ

  exit(EXIT_SUCCESS);
}