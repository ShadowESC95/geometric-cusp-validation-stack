// Geant4 integration template: map a tabulated emergent m_eff(x) field into
// a custom particle process. This is intentionally a hook, because the actual
// particle species and equation of motion must come from the physical model.
#include "G4ThreeVector.hh"
#include <vector>
#include <algorithm>
#include <cmath>

class GeometricMassField {
public:
  explicit GeometricMassField(std::vector<double> x, std::vector<double> m)
    : xs_(std::move(x)), ms_(std::move(m)) {}
  double mass(double x) const {
    if(xs_.empty()) return 0.0;
    if(x<=xs_.front()) return ms_.front();
    if(x>=xs_.back()) return ms_.back();
    auto it=std::upper_bound(xs_.begin(),xs_.end(),x);
    std::size_t j=std::size_t(it-xs_.begin()); double t=(x-xs_[j-1])/(xs_[j]-xs_[j-1]);
    return (1-t)*ms_[j-1]+t*ms_[j];
  }
private:
  std::vector<double> xs_,ms_;
};
