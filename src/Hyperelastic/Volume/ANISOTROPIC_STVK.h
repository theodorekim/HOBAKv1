/*
This file is part of HOBAK.

Permission is hereby granted to use this software solely for non-commercial applications
and purposes including academic or industrial research, evaluation and not-for-profit media
production.
 
THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PIXAR OR ITS AFFILIATES, YALE, OR
THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef VOLUME_ANISOTROPIC_STVK_H
#define VOLUME_ANISOTROPIC_STVK_H

#include "HYPERELASTIC.h"

namespace HOBAK {
namespace VOLUME {

class ANISOTROPIC_STVK : public HYPERELASTIC
{
public:
    ANISOTROPIC_STVK(const REAL& mu, const VECTOR3& a);
    ~ANISOTROPIC_STVK() {};

    // get the strain energy
    virtual REAL psi(const MATRIX3& F) const override;

    virtual MATRIX3 PK1(const MATRIX3& F) const override;

    virtual std::string name() const override;

    virtual MATRIX9 hessian(const MATRIX3& F) const override;

    virtual MATRIX9 clampedHessian(const MATRIX3& F) const override;
    
    virtual bool energyNeedsSVD() const override;

    virtual bool PK1NeedsSVD() const override;

private:
    REAL _mu;
    
    VECTOR3 _a;
    MATRIX3 _A;
};

} // VOLUME
} // HOBAK

#endif
