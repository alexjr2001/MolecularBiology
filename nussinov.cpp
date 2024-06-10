#include<iostream>
#include<vector>
#include<algorithm>
#include<limits.h>

using namespace std;

class Cell{
    public:
    int value;
    pair<int,int> prev;
    Cell(){
        value=0;
        prev={0,0};
    }
    Cell(int _v,pair<int,int> _prev={0,0}){
        value=_v;
        prev=_prev;
    }
};

int alpha(char a, char b){
     return ((a == 'G' && b == 'C') || (a == 'C' && b == 'G') ||
           (a == 'A' && b == 'U') || (a == 'U' && b == 'A'))*-1;
}

int alpha2(char a, char b){
    int ans=0;
    if ((a == 'G' && b == 'C') || (a == 'C' && b == 'G')) ans=-5;
    else if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) ans=-4;
    else if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) ans=-1;
    return ans;
}

void calculate_min(vector<vector<Cell>>& M, int i , int j, string seq, int(*alpha)(char,char)){
    pair<int,int> cur_prev={0,0};
    M[i][j].value=min({M[i+1][j].value,M[i][j-1].value,M[i+1][j-1].value+alpha(seq[i],seq[j])});
    for (int k = i+1; k < j; ++k) {
        M[i][j] = min(M[i][j].value, M[i][k].value + M[k+1][j].value);
    }

    if (M[i][j].value == M[i+1][j].value) M[i][j].prev = {i+1, j};
    else if (M[i][j].value == M[i][j-1].value) M[i][j].prev = {i, j-1};
    else if (M[i][j].value == M[i+1][j-1].value + alpha(seq[i], seq[j])) M[i][j].prev = {i+1, j-1};
    for (int k = i + 1; k < j; ++k) {
        int newVal = M[i][k].value + M[k+1][j].value;
        if (newVal == M[i][j].value) M[i][j].prev = {k, k+1};
    }
}

void fillMatrix(vector<vector<Cell>>& M, string seq){
    int j;
    for(int d=1;d<seq.size();d++){
        for(int i=0;i+d<seq.size();i++){
            j=d+i;
            calculate_min(M,i,j,seq,alpha2);
        }
    }
}

void backtracking(const vector<vector<Cell>>& M, int i , int j, string seq){
    
    cout<<"Seq "<<M[i][j].prev.first <<" "<< M[i][j].prev.second<<endl;
    if(M[i][j].prev.first == 0 && M[i][j].prev.second == 0){
        return;
    }

    cout<<endl<<"| "<<" "<<" |"<<endl;
    
    if(M[i][j].prev.first == i+1 && M[i][j].prev.second == j-1) cout<<seq[i+1]<<" - "<<seq[j-1];
    else if(M[i][j].prev.first == i+1) cout<<seq[i+1]<<"   |";
    else if(M[i][j].prev.second == j-1) cout<<"|   "<<seq[j-1];


    backtracking(M, M[i][j].prev.first, M[i][j].prev.second, seq);
}

void backtracking2(const vector<vector<Cell>>& M, int i , int j, string seq, int from){
    
    //cout<<"Seq "<<M[i][j].prev.first <<" "<< M[i][j].prev.second<<endl;
    if(M[i][j].prev.first == 0 && M[i][j].prev.second == 0){
        return;
    }

    if(from==1) cout<<seq[i]<<" - "<<seq[j];
    else if(from==2) cout<<seq[i]<<"   |";
    else if(from==0) cout<<"|   "<<seq[j];
    else cout<<seq[i]<<"   "<<seq[j];
    cout<<endl<<"| "<<" "<<" |"<<endl;

    
    if(M[i][j].prev.first == i+1 && M[i][j].prev.second == j-1) from=1;
    else if(M[i][j].prev.first == i+1) from=2;
    else if(M[i][j].prev.second == j-1) from=0;
    else from=3;


    backtracking2(M, M[i][j].prev.first, M[i][j].prev.second, seq, from);
}

void iterative_build(const vector<vector<Cell>>& M, int i, int j, string seq){
    int n = seq.size();
    i = M[0][n-1].prev.first;
    j = M[0][n-1].prev.second;
    int from = 1;

    while(i>0 || j>0){
        if(from==1) cout<<seq[i]<<" - "<<seq[j];
        else if(from==2) cout<<seq[i]<<"   |";
        else if(from==0) cout<<"|   "<<seq[j];
        else cout<<seq[i]<<"   "<<seq[j];
        cout<<endl<<"| "<<" "<<" |"<<endl;

    
        if(M[i][j].prev.first == i+1 && M[i][j].prev.second == j-1) from=1;
        else if(M[i][j].prev.first == i+1) from=2;
        else if(M[i][j].prev.second == j-1) from=0;
        else{
            from=3;
            iterative_build(M,i,M[i][j].prev.first,seq);
            iterative_build(M,M[i][j].prev.second,j,seq);
        }
    }
}

void build_sequences(const vector<vector<Cell>>& M, int i , int j, string seq){
    //cout<<seq[0]<<" - "<<seq[seq.size()-1];
    backtracking2(M,i,j,seq,1);
}

void printMatrix(const vector<vector<Cell>>& M){
    for(auto vec:M){
        for(auto x:vec){
            cout<<x.value<<" ";
        }
        cout<<endl;
    }
}


int main(){
    string sequence = "GGAAAUCC";
    sequence = "ACUCGAUUCCGAG";
    vector<vector<Cell>> matrix(sequence.size(),vector<Cell>(sequence.size(),Cell()));
    fillMatrix(matrix,sequence);
    printMatrix(matrix);
    
    build_sequences(matrix,0,sequence.size()-1,sequence);


    return 0;
}