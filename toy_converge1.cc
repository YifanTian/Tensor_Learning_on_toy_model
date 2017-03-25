#include "ITensor/all.h"
#include "fstream"
#include "string"


using namespace itensor;
using namespace std;
using std::vector;

/*
Real stat(ITensor Bl, ITensor dL, ITensor ts1, ITensor ts2,
Real *pa1, Real pa2, int train_set, Index s1, Index s2) {

  int i;
  for(int i = train_set; i<train_set+100; i++) {
    ITensor dL, ts1, ts2;

    dL.fill(0.0);
    dL.set(l(pa1[i][4]),1);

    ts1.set(s1(1),pa2[i][0]);
    ts1.set(s1(2),pa2[i][1]);
    ts2.set(s2(1),pa2[i][2]);
    ts2.set(s2(2),pa2[i][3]);
    auto fl = Bl*ts1*ts2;


    dL.fill(0.0);
    dL.set(l(pa2[i][4]),1);

    ts1.set(s1(1),pa1[i][0]);
    ts1.set(s1(2),pa1[i][1]);
    ts2.set(s2(1),pa1[i][2]);
    ts2.set(s2(2),pa1[i][3]);
    fl = Bl*ts1*ts2;
  }

  return 0;
}
*/

int main()
  {

  //Real pa1[784][2];
  Real pa1[300][7],pa2[300][7];
  //Real set1[300][3],set2[300][5];
  int i=0,train_set = 200, set_N = 300;
  Real pi = 3.1415926;
  ifstream file1,file2;
  //char filename[512] = {'s'};
  //string filename ＝ 's.txt';

  //string filename1 = "gaussion1.txt";
  //string filename2 = "gaussion2.txt";
  string filename1 = "moon1.txt";
  string filename2 = "moon2.txt";
  file1.open(filename1,ios::in);
  file2.open(filename2,ios::in);

  ofstream outfile_Bl1,outfile_Bl2,outfile_dots1,outfile_dots2;
  outfile_Bl1.open("Bl1.txt"),outfile_Bl2.open("Bl2.txt"),outfile_dots1.open("classification1.txt"),outfile_dots2.open("classification2.txt");

    i = 0;
    if(file1.fail())
      {
        cout<<"file not found."<<endl;
        file1.close();
        cin.get();
        cin.get();
      }
    else
    {
      while(!file1.eof()&&i<set_N)   //读取数据到数组,file.eof()判断文件是否为空
        {
          //if (i%2 == 0) file>>pa1[i][0];
          //else file>>pa1[i][1];
          file1>>pa1[i][0];
          file1>>pa1[i][1];
          file1>>pa1[i][2];
          file1>>pa1[i][3];
          file1>>pa1[i][4];
          file1>>pa1[i][5];
          file1>>pa1[i][6];
          i++;
          //cout<<i<<" "<<pa1[i]<<endl;
        }
    }
      file1.close(); //close file/

      i = 0;
      if(file2.fail())
        {
          cout<<"file not found."<<endl;
          file2.close();
          cin.get();
          cin.get();
        }
      else
      {
        while(!file2.eof()&&i<set_N)   //读取数据到数组,file.eof()判断文件是否为空
          {
            file2>>pa2[i][0];
            file2>>pa2[i][1];
            file2>>pa2[i][2];
            file2>>pa2[i][3];
            file2>>pa2[i][4];
            file2>>pa2[i][5];
            file2>>pa2[i][6];
            i++;
            //cout<<i<<" "<<pa1[i]<<endl;
          }
      }
        file2.close(); //close file/

  //for(i=0;i<783;i++)
  for(int j=0;j<100;j++)
  {
    cout<<j<<" "<<pa1[j][0]<<" "<<pa1[j][1]<<" "<<pa1[j][2]<<" "<<pa1[j][3]<<" "<<pa1[j][4]<<" "<<pa1[j][5]<<" "<<pa1[j][6]<<endl;
    cout<<j<<" "<<pa2[j][0]<<" "<<pa2[j][1]<<" "<<pa2[j][2]<<" "<<pa2[j][3]<<" "<<pa2[j][4]<<" "<<pa2[j][5]<<" "<<pa2[j][6]<<endl;
  }

  for(int i = 1; i<20; i++) {
    cout<<pa2[i][5]<<" "<<pa2[i][6]<<" "<<pa2[i][4]<<" "<<endl;
  }

  int N = 2;
  auto sites = SpinHalf(N);
  auto state = InitState(sites);          //??
  auto psi = MPS(state);

  auto l = Index("classification",2);
  auto s1 = Index("x1",2);
  auto s2= Index("y2",2);

  auto Bl = ITensor(s1,s2,l);
  randomize(Bl);
  auto ts1 = ITensor(s1);
  auto ts2 = ITensor(s2);
  auto ts = ITensor(s1,s2);
  auto dL = ITensor(l);
  auto fl = ITensor(l);
  auto dBl = ITensor(s1,s2,l);

  auto cost = 0.0;
  for(int j = 1; j<20; j++) {
    dL.fill(0.0);
    dL.set(l(pa1[j][4]),1);
    ts1.set(s1(1),pa2[j][0]); ts1.set(s1(2),pa2[j][1]);
    ts2.set(s2(1),pa2[j][2]); ts2.set(s2(2),pa2[j][3]);
    auto fl = Bl*ts1*ts2;
    cost += pow(norm(dL-fl),2);
    //PrintData(cost);

    dL.fill(0.0);
    dL.set(l(pa2[j][4]),1);
    ts1.set(s1(1),pa1[j][0]); ts1.set(s1(2),pa1[j][1]);
    ts2.set(s2(1),pa1[j][2]); ts2.set(s2(2),pa1[j][3]);
    fl = Bl*ts1*ts2;
    cost += pow(norm(dL-fl),2);
    //PrintData(cost);
  }
  PrintData(cost);
  cost = 0.0;

  for (int bl1 = 1;bl1 <= 2;bl1++) {
    for (int bl2 = 1;bl2 <= 2;bl2++) {
      for (int bl3 = 1;bl3 <= 2;bl3++) {
        //cout<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<endl;
        outfile_Bl1<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";
      }
    }
  }
  outfile_Bl1<<endl;

  auto outcome = 0;
  for(int i = 1; i<20; i++) {

  ts1.set(s1(1),pa2[i][0]);
  ts1.set(s1(2),pa2[i][1]);
  ts2.set(s2(1),pa2[i][2]);
  ts2.set(s2(2),pa2[i][3]);

  fl = Bl*ts1*ts2;
  outfile_dots1<<pa2[i][5]<<" "<<pa2[i][6]<<" "<<pa2[i][4]<<" ";
  auto result = 0;
  if (fl.real(l(1))>fl.real(l(2))) {
    result = 1;
    outfile_dots1<<result; }
  else {
    result = 2;
    outfile_dots1<<result; }
  if (result == pa2[i][4]) outcome += 1;
  outfile_dots1<<endl;
  //PrintData(fl);

  ts1.set(s1(1),pa1[i][0]);
  ts1.set(s1(2),pa1[i][1]);
  ts2.set(s2(1),pa1[i][2]);
  ts2.set(s2(2),pa1[i][3]);

  fl = Bl*ts1*ts2;
  //cout<<"i 2 "<<i<<endl;
  //cout<<i<<" "<<pa1[i][0]<<" "<<pa1[i][1]<<" "<<pa1[i][2]<<" "<<pa1[i][3]<<" "<<pa1[i][4]<<endl;
  //cout<<i<<" "<<pa1[i][5]<<" "<<pa1[i][6]<<endl;
  outfile_dots1<<pa1[i][5]<<" "<<pa1[i][6]<<" "<<pa1[i][4]<<" ";
  result = 0;
  if (fl.real(l(1))>fl.real(l(2))) {
    result = 1;
    outfile_dots1<<result; }
  else {
    result = 2;
    outfile_dots1<<result; }
  if (result == pa1[i][4]) outcome += 1;
  outfile_dots1<<endl;
}
outfile_dots1<<outcome;
cout<<"original classification: "<<outcome<<endl;

for(int i = 1; i<20; i++) {
  cout<<pa2[i][5]<<" "<<pa2[i][6]<<" "<<pa2[i][4]<<" "<<endl;
}


  //for(int i = 0; i<train_set; i++) {
for (int sweep = 1; sweep < 50; sweep++)  {

  dBl.fill(0.0);
  for(int i = 0; i<100; i++) {
    dL.fill(0.0);
    dL.set(l(pa1[i][4]),1);

    ts1.set(s1(1),pa1[i][0]);
    ts1.set(s1(2),pa1[i][1]);
    ts2.set(s2(1),pa1[i][2]);
    ts2.set(s2(2),pa1[i][3]);

    fl = Bl*ts1*ts2;
    dBl = dBl + 0.01*(ts1*ts2*(dL-fl));              //gradient

    //cost = stat(Bl,dL,ts1,ts2,&pa1,pa2,train_set,s1,s2);
    int j;
    for(int j = train_set; j<train_set+50; j++) {
      dL.fill(0.0);
      dL.set(l(pa1[j][4]),1);
      ts1.set(s1(1),pa2[j][0]); ts1.set(s1(2),pa2[j][1]);
      ts2.set(s2(1),pa2[j][2]); ts2.set(s2(2),pa2[j][3]);
      auto fl = Bl*ts1*ts2;
      cost += pow(norm(dL-fl),2);
      //PrintData(cost);

      dL.fill(0.0);
      dL.set(l(pa2[j][4]),1);
      ts1.set(s1(1),pa1[j][0]); ts1.set(s1(2),pa1[j][1]);
      ts2.set(s2(1),pa1[j][2]); ts2.set(s2(2),pa1[j][3]);
      fl = Bl*ts1*ts2;
      cost += pow(norm(dL-fl),2);
      //PrintData(cost);
    }
    //cout<<i<<" "<<"cost = "<<cost<<endl;
    cost = 0.0;

    /*
    int bln = 0;
    for (int bl1 = 1;bl1 <= 2;bl1++) {
      for (int bl2 = 1;bl2 <= 2;bl2++) {
        for (int bl3 = 1;bl3 <= 2;bl3++) {
          bln++;
          //cout<<bl1<<" "<<bl2<<" "<<bl3<<" "<<bln<<" "<<4*(bl1-1)+2*(bl2-1)+(bl3-1)+1<<" ";
          //cout<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<endl;
          outfile_Bl2<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";
        }
      }
    }
    outfile_Bl2<<endl;
    */

    dL.fill(0.0);
    dL.set(l(pa2[i][4]),1);

    ts1.set(s1(1),pa2[i][0]);
    ts1.set(s1(2),pa2[i][1]);
    ts2.set(s2(1),pa2[i][2]);
    ts2.set(s2(2),pa2[i][3]);

    ts = ts1*ts2;
    fl = Bl*ts1*ts2;
    dBl = dBl + 0.01*(ts1*ts2*(dL-fl));              //gradient
    //Bl+=dBl;
                                            //calculate cost function
    for(int j = train_set; j<train_set+50; j++) {
      dL.fill(0.0);
      dL.set(l(pa1[j][4]),1);
      ts1.set(s1(1),pa2[j][0]); ts1.set(s1(2),pa2[j][1]);
      ts2.set(s2(1),pa2[j][2]); ts2.set(s2(2),pa2[j][3]);
      auto fl = Bl*ts1*ts2;
      cost += pow(norm(dL-fl),2);
      //PrintData(cost);

      dL.fill(0.0);
      dL.set(l(pa2[j][4]),1);
      ts1.set(s1(1),pa1[j][0]); ts1.set(s1(2),pa1[j][1]);
      ts2.set(s2(1),pa1[j][2]); ts2.set(s2(2),pa1[j][3]);
      fl = Bl*ts1*ts2;
      cost += pow(norm(dL-fl),2);
      //PrintData(cost);
    }
    //cout<<i<<" "<<"cost = "<<cost<<endl;
    cost = 0.0;

    /*
      for (int bl1 = 1;bl1 <= 2;bl1++) {
        for (int bl2 = 1;bl2 <= 2;bl2++) {
          for (int bl3 = 1;bl3 <= 2;bl3++) {
            //cout<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<endl;
            outfile_Bl2<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";
          }
        }
      }
      outfile_Bl2<<endl;
    */

  }

  cost = 0.0;
  Bl = Bl + dBl;
  /*
  for(int j = 1; j<20; j++) {
    dL.fill(0.0);
    dL.set(l(pa1[j][4]),1);
    ts1.set(s1(1),pa2[j][0]); ts1.set(s1(2),pa2[j][1]);
    ts2.set(s2(1),pa2[j][2]); ts2.set(s2(2),pa2[j][3]);
    auto fl = Bl*ts1*ts2;
    cost += pow(norm(dL-fl),2);
    //PrintData(cost);

    dL.fill(0.0);
    dL.set(l(pa2[j][4]),1);
    ts1.set(s1(1),pa1[j][0]); ts1.set(s1(2),pa1[j][1]);
    ts2.set(s2(1),pa1[j][2]); ts2.set(s2(2),pa1[j][3]);
    fl = Bl*ts1*ts2;
    cost += pow(norm(dL-fl),2);
    //PrintData(cost);
  }
  //PrintData(cost);
  cout<<sweep<<" "<<"cost = "<<cost<<endl;
  */
}

for (int bl1 = 1;bl1 <= 2;bl1++) {
  for (int bl2 = 1;bl2 <= 2;bl2++) {
    for (int bl3 = 1;bl3 <= 2;bl3++) {
      //cout<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<endl;
      outfile_Bl2<<Bl.real(s1(bl1),s2(bl2),l(bl3))<<" ";
    }
  }
}
outfile_Bl2<<endl;

cost = 0.0;
for(int j = 1; j<20; j++) {
  dL.fill(0.0);
  dL.set(l(pa1[j][4]),1);
  ts1.set(s1(1),pa2[j][0]); ts1.set(s1(2),pa2[j][1]);
  ts2.set(s2(1),pa2[j][2]); ts2.set(s2(2),pa2[j][3]);
  auto fl = Bl*ts1*ts2;
  cost += pow(norm(dL-fl),2);
  //PrintData(cost);

  dL.fill(0.0);
  dL.set(l(pa2[j][4]),1);
  ts1.set(s1(1),pa1[j][0]); ts1.set(s1(2),pa1[j][1]);
  ts2.set(s2(1),pa1[j][2]); ts2.set(s2(2),pa1[j][3]);
  fl = Bl*ts1*ts2;
  cost += pow(norm(dL-fl),2);
  //PrintData(cost);
}
PrintData(cost);

  outcome = 0;
  //for(int i = train_set; i<set_N; i++) {
  //for(int i = train_set; i<train_set+100; i++) {
    for(int i = 1; i<20; i++) {

    ts1.set(s1(1),pa2[i][0]);
    ts1.set(s1(2),pa2[i][1]);
    ts2.set(s2(1),pa2[i][2]);
    ts2.set(s2(2),pa2[i][3]);

    fl = Bl*ts1*ts2;
    cout<<pa2[i][5]<<" "<<pa2[i][6]<<" "<<pa2[i][4]<<" "<<endl;
    outfile_dots2<<pa2[i][5]<<" "<<pa2[i][6]<<" "<<pa2[i][4]<<" ";
    auto result = 0;
    if (fl.real(l(1))>fl.real(l(2))) {
      result = 1;
      outfile_dots2<<result; }
    else {
      result = 2;
      outfile_dots2<<result; }
    if (result == pa2[i][4]) outcome += 1;
    outfile_dots2<<endl;
    //PrintData(fl);

    ts1.set(s1(1),pa1[i][0]);
    ts1.set(s1(2),pa1[i][1]);
    ts2.set(s2(1),pa1[i][2]);
    ts2.set(s2(2),pa1[i][3]);

    fl = Bl*ts1*ts2;
    //cout<<"i 2 "<<i<<endl;
    //cout<<i<<" "<<pa1[i][0]<<" "<<pa1[i][1]<<" "<<pa1[i][2]<<" "<<pa1[i][3]<<" "<<pa1[i][4]<<endl;
    //cout<<i<<" "<<pa1[i][5]<<" "<<pa1[i][6]<<endl;
    outfile_dots2<<pa1[i][5]<<" "<<pa1[i][6]<<" "<<pa1[i][4]<<" ";
    result = 0;
    if (fl.real(l(1))>fl.real(l(2))) {
      result = 1;
      outfile_dots2<<result; }
    else {
      result = 2;
      outfile_dots2<<result; }
    if (result == pa1[i][4]) outcome += 1;
    outfile_dots2<<endl;
  }
  outfile_dots2<<outcome;
  cout<<"new classification: "<<outcome<<endl;

  outfile_Bl1.close();
  outfile_Bl2.close();

  outfile_dots1.close();
  outfile_dots2.close();

    return 0;
  }
