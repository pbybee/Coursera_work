#include <iostream>
#include <vector>
#include <algorithm>

using std::vector;
using std::cin;
using std::cout;
using std::swap;
using std::pair;
using std::make_pair;

class HeapBuilder {
 private:
  vector<int> data_;
  vector< pair<int, int> > swaps_;
  int largest;

  void WriteResponse() const {
    cout << swaps_.size() << "\n";
    for (int i = 0; i < swaps_.size(); ++i) {
      cout << swaps_[i].first << " " << swaps_[i].second << "\n";
    }
  }

  void ReadData() {
    int n;
    cin >> n;
    data_.resize(n);
    for(int i = 0; i < n; ++i)
      cin >> data_[i];
  }

  void GenerateSwaps() {
    swaps_.clear();
    
    for (int i = data_.size()/2; i >= 1; i--){
        SiftDown(i);
    }
    
//    for (int i = 0; i < data_.size(); ++i)
//      for (int j = i + 1; j < data_.size(); ++j) {
//        if (data_[i] > data_[j]) {
//          swap(data_[i], data_[j]);
//          swaps_.push_back(make_pair(i, j));
//        }
  }
  
void SiftDown(int i){
    int maxIdx = i;
    int l = Left(i);
    if (l <= data_.size() && data_[l] > data_[maxIdx]){
        maxIdx = l;
    }
    int r = Right(i);
    if (r <= data_.size() && data_[r] > data_[maxIdx]){
        maxIdx = r;
    }
    if (i != maxIdx){
        swaps_.push_back(make_pair(i,maxIdx));
        SiftDown(maxIdx);
    }
}

int Left(int i){
    return(2*i);
}

int Right(int i){
    return(2*i + 1);
}

 public:
  void Solve() {
    ReadData();
    GenerateSwaps();
    WriteResponse();
  }
};

int main() {
  std::ios_base::sync_with_stdio(false);
  HeapBuilder heap_builder;
  heap_builder.Solve();
  return 0;
}
