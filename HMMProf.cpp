#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <fstream>

using namespace std;


//Auxiliar functions
bool is_conserve(int idx, set<int> no_col){          //If is conserve region the column
    return no_col.find(idx)==no_col.end()? true:false;
}

void handle_conservative_regions(vector<string> seq, vector<string>& HMM, int i, int& j, set<int> no_col){      //Put in HMM what type of edges it has (Match,Delete,ComebackDelete(T))
    if (seq[i][j-1]!='-'&&seq[i][j]!='-') HMM[i].push_back('M');            //If match
    else if(seq[i][j-1]=='-'&&seq[i][j]!='-') HMM[i].push_back('T');        //If came back from deletion
    else if(seq[i][j-1]!='-'&&seq[i][j]=='-'){
        while(seq[i][j]=='-' && is_conserve(j,no_col)){             //Delete until a letter as long as is conservative region
            HMM[i].push_back('D');
            j++;
        }
        j--;
    } 
}

void print_legend(){
    cout<<"Our output is divided by levels, where explains the probability of every edge, here's the legend:\n\n";
    cout<<"M->Match (it goes to the next state in the main branch)\n";
    cout<<"D->Delete (New delete state out of branch)\n";
    cout<<"T->DeleteRejoin (Delete state rejoins the main branch)\n";
    cout<<"I->Insert (Creates a new node of inserting)\n";
    cout<<"J->No-Insert (what is not inserting goes straight to next state in main)\n";
    cout<<"L->Loop (Loop in the insert state)\n";
    cout<<"R->Rejoin Insertion (From insertion node returns to main branch)\n\n\n";

    return;
}


vector<string> constructHMM(vector<string> seq, set<int> no_col){
    vector<string> HMM(seq.size(),"");  //Match->M , Delete->D , Insert->I , Loop -> L , 
                                        //RejoinInsertion->R, RejoinDelete -> T, Jump -> J
    char pre_non_conservative;
    for (int i = 0; i < seq.size(); i++)
    {
        for (int j = 1; j < seq[i].size(); j++)
        {
            if(is_conserve(j-1,no_col) && is_conserve(j,no_col)){  //If both are conservative region
                handle_conservative_regions(seq,HMM,i,j,no_col);
            }
            else if(is_conserve(j-1,no_col) && !is_conserve(j,no_col)){ //Begin non-conservative(If is not conservative region, but previous is)
                pre_non_conservative=HMM[i][HMM[i].size()-1];  //Last letter pre-non-conserve region
                if (seq[i][j]=='-'){
                    while(!is_conserve(j,no_col)){      //Loop until conserve region JUMP
                        HMM[i].push_back('J');
                        j++;
                    }
                    j--;
                }
                else if (seq[i][j]!='-')  HMM[i].push_back('I');
            }
            else if(!is_conserve(j-1,no_col) && !is_conserve(j,no_col)){    //Previos and current are not conservatives regions
                if(seq[i][j]=='-'){ 
                    while(!is_conserve(j,no_col)){     
                        HMM[i].push_back('R');
                        j++;
                    }
                    j--;
                }
                else if(seq[i][j]!='-') HMM[i].push_back('L');
            }
            else if(!is_conserve(j-1,no_col) && is_conserve(j,no_col)){ //If is not conservative region, but post is
                if(seq[i][j]=='-'){
                    while(seq[i][j]=='-' && is_conserve(j,no_col)){
                        HMM[i].push_back('D');
                        j++;
                    }
                    j--;
                } 
                else{
                    if (pre_non_conservative=='D') HMM[i].push_back('T');       //If there's a letter an last letter prev-conserve is delete so we join it
                    else HMM[i].push_back('M');
                } 
            }
        }
    }

    return HMM;
}

void printHMM(vector<string> HMM){
    for (auto i:HMM){
        cout<<i<<endl;
    }
}

set<int> conservativeRegion(vector<string> seq){        //Which columns have more than the half residuals and not hyphen
    int count=0;
    set<int> no_columns;
    for(int j=0;j<seq[0].size();j++){
        for(auto i:seq){
            if(i[j]!='-') count++;
        }
        if(count<=seq.size()/2){
            no_columns.emplace(j);
        }
        count=0;
    }
    return no_columns;
}

void calculateProb(vector<string> HMM, set<int> no_col){
    double total=HMM.size();
    map<char,double> times;
    int level=0;
    for(int j=0; j<HMM[0].size();j++){
        for (auto i:HMM)            //Count how many occurrences in the level
        {
            if(times[i[j]]) times[i[j]]++;
            else times[i[j]]=1;
        }

        //If we get into a no-conserve region so we print M and I, but we count Jumps Insertions, Returns.
        //And at the end we print which insertion were in loop and how much came back
        if(is_conserve(j,no_col)){      //In case we are in conserve region print D or M
            cout<<"Level "<<level<<endl;
            for(auto [key,value]:times){    
                cout<<key<<": "<<value/total<<endl;
            }
            if(j+1<HMM[0].size() && is_conserve(j+1,no_col)) times.clear();
            level++;
        }
        else if(j+1<HMM[0].size() && is_conserve(j+1,no_col)){  //Back to conserve region
            double no_conserve_total=times['L']+times['R']+times['I'];
            cout<<"L: "<< times['L']/no_conserve_total<<endl;       //Times in loop
            cout<<"R: "<< 1-times['L']/no_conserve_total<<endl;  //Times of insertion or rejoin
            times.clear();
        }
    }
}

int main() {

    //Init sequence alignment
    vector<string> seq{
        "VGA--HAGEY",
        "V----NVDEV",
        "VEA--DVAGH",
        "VKG------D",
        "VYS--TYETS",
        "FNA--NIPKH",
        "IAGADNGAGY"
    };
    seq={
        "ACA---ATG",
        "TCAACTATC",
        "ACAC--AGC",
        "ACA---ATC",
        "A-C---ATC"
    };

    //Select which columns are not conservatives
    set<int> no_columns=conservativeRegion(seq);

    //Construct which edges has every column in conservative region
    vector<string> HMM=constructHMM(seq,no_columns);
    //printHMM(HMM);

    print_legend();
    calculateProb(HMM,no_columns);
    

    return 0;
}