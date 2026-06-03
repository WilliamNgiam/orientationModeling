#include <module/Module.h>
#include <distribution/ScalarDist.h>

#include "distributions/DVonMises.h"

namespace jags {
namespace vonmises {

class VONMISESModule : public Module {
 public:
  VONMISESModule();
  ~VONMISESModule();
};

VONMISESModule::VONMISESModule() : Module("vonmises") {
  insert(new DVonMises);
}

VONMISESModule::~VONMISESModule() {
  std::vector<Distribution *> const &dvec = distributions();
  for (unsigned int i = 0; i < dvec.size(); ++i) {
    delete dvec[i];
  }
}

}  // namespace vonmises
}  // namespace jags

jags::vonmises::VONMISESModule _vonmises_module;
