#include <psi4/libmints/molecule.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libmints/wavefunction.h>
#include <psi4/libmints/factory.h>
#include <psi4/libmints/integral.h>
#include <psi4/libmints/electrostatic.h>
#include <psi4/libmints/oeprop.h>
#include <psi4/libpsi4util/process.h>


using namespace psi;

std::vector<double> calculate_esp_at_points(std::shared_ptr<Wavefunction> wfn, std::vector<Vector3> points) {
    std::shared_ptr<BasisSet> basisset = wfn->basisset();
    std::shared_ptr<Molecule> mol = basisset->molecule();
    std::shared_ptr<IntegralFactory> integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset, basisset, basisset, basisset));
    std::shared_ptr<ElectrostaticInt> epot(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));
    auto oeprop(std::make_shared<OEProp>(wfn));

    int n_atoms = mol->natom();
    int nbf = basisset->nbf();
    std::vector<double> esp_values;

    SharedMatrix Dtot = oeprop->Da_ao();
    if (wfn->same_a_b_dens()) {
        Dtot->scale(2.0);
    } else {
        Dtot->add(oeprop->Db_ao());
    }

    // compute the electrostatic potential at each of the points
    // this code probably should be in a different function?
    outfile->Printf("\n Electrostatic potentials at van der Waals shells:\n");
    outfile->Printf(" ---------------------------------------------------\n");
    outfile->Printf("   x   y   z  Electrostatic Potential (a.u.)\n");
    outfile->Printf(" ---------------------------------------------------\n");
    for (size_t i = 0; i < points.size(); i++) {
        std::stringstream s;
        SharedMatrix ints(new Matrix(s.str(), nbf, nbf));
        epot->compute(ints, points[i]);
        double elec = Dtot->vector_dot(ints);
        double nuc = 0.0;

        for (int atom1 = 0; atom1 < n_atoms; atom1++)
            nuc += (mol->Z(atom1) / points[i].distance(mol->xyz(atom1)));

        esp_values.push_back(nuc+elec);
        outfile->Printf("     %8.5f %8.5f %8.5f    %16.12f\n", points[i][0], points[i][1], points[i][2], nuc+elec);
    }
    outfile->Printf(" ---------------------------------------------------\n");

    return esp_values;
}
