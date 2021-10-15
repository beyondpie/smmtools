#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// From ArchR General_Utils.cpp
//' @export
// [[Rcpp::export]]
IntegerMatrix tabulate2dCpp(IntegerVector &x, int xmin, int xmax, IntegerVector &y, int ymin, int ymax){
  if(x.size() != y.size()){
    stop("width must equal size!");
  }
  int n = x.size();
  IntegerVector rx = seq(xmin, xmax);
  IntegerVector ry = seq(ymin, ymax);
  IntegerMatrix mat( ry.size() , rx.size() );
  int rys = ry.size();
  int rxs = rx.size();
  int xi,yi;
  for(int i = 0; i < n; i++){
    xi = (x[i] - xmin);
    yi = (y[i] - ymin);
    if(yi >= 0 && yi < rys){
      if(xi >= 0 && xi < rxs){
        mat( yi , xi ) = mat( yi , xi ) + 1; 
      }
    }
  }
  return mat;
}
