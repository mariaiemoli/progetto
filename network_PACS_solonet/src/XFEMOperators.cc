#include "../include/XFEMOperators.h"

namespace getfem
{

// Defining unit normal on a level set ------------------------------------

level_set_unit_normal::level_set_unit_normal ( const getfem::mesh_fem& mf_,
                                               const scalarVector_Type& U_ ) :
    mf(mf_), U(mf_.nb_basic_dof()), N(mf_.linked_mesh().dim()), gradU(1, N)
{
    sizes_.resize(1);
    sizes_ [ 0 ] = short_type(N);
    mf.extend_vector(U_, U);
}

void level_set_unit_normal::compute ( getfem::fem_interpolation_context& ctx,
                                      bgeot::base_tensor& t )
{
    size_type cv = ctx.convex_num();
    coeff.resize(mf.nb_basic_dof_of_element(cv));
    gmm::copy(gmm::sub_vector(U,
            gmm::sub_index(mf.ind_basic_dof_of_element(cv))), coeff);
    ctx.pf()->interpolation_grad(ctx, coeff, gradU, 1);
    scalar_type norm = gmm::vect_norm2(gmm::mat_row(gradU, 0));
    for ( size_type i = 0; i < N; ++i )
    {
        t [ i ] = gradU(0, i) / norm;
    }
}

void interpolateVelocity ( scalarVector_Type& V,
                           const base_node& point,
                           const size_type &elem,
                           const scalarVector_Type& gdl,
                           const MeshHandlerPtr_Type& mesh,
                           const sparseMatrixPtr_Type& globalMatrix )
{
    // componenti x,y, punti in cui voglio v (x,y), elementi interessati,
    // matrice del sistema per dedurre i segni, shift per i dl raddoppiati,
    // vettore dei gradi di libertà

    scalar_type segno, lato, coeff;
    size_type idof;
    const size_type shiftVelocity = gdl.size();

    // area dell'elemento corrispondente
    const scalar_type area = mesh->getMesh().convex_area_estimate(elem);

    // vertici dell'elemento
    base_node xx1 = mesh->getMesh().points_of_convex(elem) [ 0 ];
    base_node xx2 = mesh->getMesh().points_of_convex(elem) [ 1 ];
    base_node xx3 = mesh->getMesh().points_of_convex(elem) [ 2 ];

    //prima funzione di base
    idof = mesh->getMeshFEMVector().ind_basic_dof_of_element(elem) [ 0 ];

    //per dedurre il segno degil rt0 - convenzione per la normale
    segno = -(*globalMatrix)(idof,
            mesh->getMeshFEMScalar().ind_basic_dof_of_element(elem) [ 0 ]
                    + shiftVelocity);

    /*if ( segno == 0 ) //se il segno è zero vuol dire che sto guardando un elemento tagliato, vado a prendere il gdl esteso di pressione
     for ( size_type j = 0; j < M_mesh->getExtendedDOFScalar().size(); ++j )
     {
     if ( M_mesh->getExtendedDOFScalar(j) == i )
     {
     segno = -(*M_globalMatrix)(idof,
     M_mesh->getMeshFEMScalar().nb_dof() + shiftVelocity + j);
     }
     }*/

    segno /= (gmm::abs(segno) > 0) ? gmm::abs(segno) : 1;
    lato = std::pow(std::pow(xx3 [ 0 ] - xx2 [ 0 ], 2) + std::pow(xx3 [ 1 ]
            - xx2 [ 1 ], 2), 0.5); //calcolo la lunghezza del lato

    //sommo il contributo della funzione di base
    coeff = gdl [ idof ] / (2 * area) * segno * lato;
    V [ 0 ] += coeff * (point [ 0 ] - xx1 [ 0 ]);
    V [ 1 ] += coeff * (point [ 1 ] - xx1 [ 1 ]);

    //seconda funzione di base
    idof = mesh->getMeshFEMVector().ind_basic_dof_of_element(elem) [ 1 ];
    segno = -(*globalMatrix)(idof,
            mesh->getMeshFEMScalar().ind_basic_dof_of_element(elem) [ 0 ]
                    + shiftVelocity);
    /*
     if ( segno == 0 )
     for ( size_type j = 0; j < M_mesh->getExtendedDOFScalar().size(); ++j )
     {
     if ( M_mesh->getExtendedDOFScalar(j) == i )
     {
     segno = -(*M_globalMatrix)(idof,
     M_mesh->getMeshFEMScalar().nb_dof() + shiftVelocity + j);
     }
     }
     */
    segno /= (gmm::abs(segno) > 0) ? gmm::abs(segno) : 1;
    lato = std::pow(std::pow(xx3 [ 0 ] - xx1 [ 0 ], 2) + std::pow(xx3 [ 1 ]
            - xx1 [ 1 ], 2), 0.5);

    coeff = gdl [ idof ] / (2 * area) * segno * lato;
    V [ 0 ] += coeff * (point [ 0 ] - xx2 [ 0 ]);
    V [ 1 ] += coeff * (point [ 1 ] - xx2 [ 1 ]);

    //terza funzione di base
    idof = mesh->getMeshFEMVector().ind_basic_dof_of_element(elem) [ 2 ];
    segno = -(*globalMatrix)(idof,
            mesh->getMeshFEMScalar().ind_basic_dof_of_element(elem) [ 0 ]
                    + shiftVelocity);

    /*if ( segno == 0 )
     {
     for ( size_type j = 0; j < M_mesh->getExtendedDOFScalar().size(); ++j )
     {
     if ( M_mesh->getExtendedDOFScalar(j) == i )
     {
     segno = -(*M_globalMatrix)(idof,
     M_mesh->getMeshFEMScalar().nb_dof() + shiftVelocity + j);
     }
     }
     }*/

    segno /= (gmm::abs(segno) > 0) ? gmm::abs(segno) : 1;
    lato = std::pow(std::pow(xx2 [ 0 ] - xx1 [ 0 ], 2) + std::pow(xx2 [ 1 ]
            - xx1 [ 1 ], 2), 0.5);

    coeff = gdl [ idof ] / (2 * area) * segno * lato;
    V [ 0 ] += coeff * (point [ 0 ] - xx3 [ 0 ]);
    V [ 1 ] += coeff * (point [ 1 ] - xx3 [ 1 ]);

} // interpolateVelocity

//calcola la lunghezza dei lati della frattura
//integrando una funzione di base p0 ottengo l'area, o la lunghezza
void stimaLati ( scalarVector_Type& latiF,
                 const FractureHandlerPtr_Type& fracture )
{
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();

    getfem::generic_assembly assem;

    assem.set("t=comp(Base(#1));"
        "V$1(#1)+=t(:);");

    assem.push_mi(fracture->getIntegrationMethodPressureVisualization());

    assem.push_mf(fracture->getMeshFEMPressureVisualization());

    assem.push_vec(latiF);

    assem.assembly(-1);

} // stima_lati

void interpolationMatrix ( sparseMatrixPtr_Type& MM,
                           const MeshHandlerPtr_Type& mesh,
                           const FractureHandlerPtr_Type& fracture,
                           const size_type& regionID )
{
    const scalar_type toll = 1.0e-9;
    base_node nodo1, nodo2;

    const scalar_type gauss_points [ ] =
    {
   1.000000000000000,
   0.982572296604548,
   0.941976296959745,
   0.879294755323590,
   0.796001926077712,
   0.694051026062223,
   0.575831960261831,
   0.444115783279002,
   0.301989856508765,
   0.152785515802185,
                   0,
  -0.152785515802185,
  -0.301989856508765,
  -0.444115783279002,
  -0.575831960261831,
  -0.694051026062223,
  -0.796001926077712,
  -0.879294755323590,
  -0.941976296959745,
  -0.982572296604548,
  -1.000000000000000
    };

    const scalar_type gauss_weights [ ] =
    {
   0.004761904761905,
   0.029184840098506,
   0.051843169000850,
   0.073273918185074,
   0.092985467957886,
   0.110517083219123,
   0.125458121190869,
   0.137458462860041,
   0.146236862447977,
   0.151587575111681,
   0.153385190332175,
   0.151587575111681,
   0.146236862447977,
   0.137458462860041,
   0.125458121190869,
   0.110517083219123,
   0.092985467957886,
   0.073273918185074,
   0.051843169000850,
   0.029184840098506,
   0.004761904761905
    };

    const size_type numberPoints = 20;

    const getfem::mesh_fem& meshFEMFracture1 =
            fracture->getMeshFEMPressureVisualization();

    const getfem::mesh_fem& meshFEMFracture2 = fracture->getMeshFEMLinear();

    const sizeVector_Type extendedDOF = mesh->getExtendedDOFScalar(
                                        fracture->getId());
    const size_type shiftExtended = extendedDOF.size();

    mesh_region regionMesh = mesh->getMesh().region(regionID);

    const size_type shiftScalar = meshFEMFracture1.nb_dof();
    for ( size_type i = 0; i < shiftScalar; ++i )
    {
        const size_type node = meshFEMFracture1.first_convex_of_basic_dof(i);

        nodo1 = meshFEMFracture2.point_of_basic_dof(node, 0);
        nodo2 = meshFEMFracture2.point_of_basic_dof(node, 1);

        for ( size_type kk = 0; kk < numberPoints; ++kk )
        {

            const base_node xx = 0.5 * (nodo1 + nodo2) + 0.5 * (nodo2 - nodo1)
                    * gauss_points [ kk ];

            for ( size_type jj = 0; jj < shiftExtended; ++jj )
            {
                const size_type j = extendedDOF [ jj ];
//                const size_type el = (regionMesh.index())[e];
                size_type el =
                                mesh->getMeshFEMCoefficients().first_convex_of_basic_dof(j);

                const base_node x0 = xx
                        - mesh->getMesh().points_of_convex(el) [ 0 ];
                const base_node x1 = xx
                        - mesh->getMesh().points_of_convex(el) [ 1 ];
                const base_node x2 = xx
                        - mesh->getMesh().points_of_convex(el) [ 2 ];

                const scalar_type A = 0.5 * (gmm::abs(x0 [ 0 ] * x1 [ 1 ]
                        - x0 [ 1 ] * x1 [ 0 ]) + gmm::abs(x0 [ 0 ] * x2 [ 1 ]
                        - x0 [ 1 ] * x2 [ 0 ]) + gmm::abs(x2 [ 0 ] * x1 [ 1 ]
                        - x2 [ 1 ] * x1 [ 0 ]));


                const scalar_type area = gmm::abs(mesh->getMesh().convex_area_estimate(el));
                if ( gmm::abs(A - area) / area <= toll )
                {
                    (*MM)(j, i) += gauss_weights [ kk ];
                    continue;
                }
            }

        }
    }

    for ( size_type j = 0; j < mesh->getMeshFEMCoefficients().nb_dof(); ++j )
    {
        scalar_type s = 0;
        for ( size_type i = 0; i < meshFEMFracture1.nb_dof(); ++i )
        {
            s += (*MM)(j, i);
        }

        if ( s > 0 )
        {
            for ( size_type i = 0; i < meshFEMFracture1.nb_dof(); ++i )
            {
                (*MM)(j, i) = (*MM)(j, i) / s;
            }
        }
    }

} // interpolationMatrix

void interpolationMatrix ( sparseMatrixPtrContainer_Type& interpolationMatrices,
                           const MeshHandlerPtr_Type& mesh,
                           FracturePtrContainer_Type& fractures,
                           const IntersectData_Type& intersect )
{
    const scalar_type toll = 1.0e-9;
    base_node nodo1, nodo2;

    const scalar_type gauss_points [ ] =
    {
   1.000000000000000,
   0.982572296604548,
   0.941976296959745,
   0.879294755323590,
   0.796001926077712,
   0.694051026062223,
   0.575831960261831,
   0.444115783279002,
   0.301989856508765,
   0.152785515802185,
                   0,
  -0.152785515802185,
  -0.301989856508765,
  -0.444115783279002,
  -0.575831960261831,
  -0.694051026062223,
  -0.796001926077712,
  -0.879294755323590,
  -0.941976296959745,
  -0.982572296604548,
  -1.000000000000000
    };

    const scalar_type gauss_weights [ ] =
    {
   0.004761904761905,
   0.029184840098506,
   0.051843169000850,
   0.073273918185074,
   0.092985467957886,
   0.110517083219123,
   0.125458121190869,
   0.137458462860041,
   0.146236862447977,
   0.151587575111681,
   0.153385190332175,
   0.151587575111681,
   0.146236862447977,
   0.137458462860041,
   0.125458121190869,
   0.110517083219123,
   0.092985467957886,
   0.073273918185074,
   0.051843169000850,
   0.029184840098506,
   0.004761904761905
    };

    const size_type numberPoints = 20;

    size_type el = intersect.getElementID();
    size_type j = mesh->getMeshFEMScalar().ind_basic_dof_of_element(el)[0];

    for ( size_type f = 0; f < fractures.size(); ++f )
    {
        for ( size_type i = 0; i < fractures[f]->getMeshFEMPressure().nb_dof() + fractures[f]->getNumExtendedPressure(); ++i )
        {
            (*(interpolationMatrices[f]))(j, i) = 0;
        }
    }

    for ( size_type f = 0; f < fractures.size(); ++f )
    {
        const getfem::mesh_fem& meshFEMFracture1 =
        fractures[f]->getMeshFEMPressureVisualization();

        const getfem::mesh_fem& meshFEMFracture2 = fractures[f]->getMeshFEMLinear();

        const size_type shiftScalar = meshFEMFracture1.nb_dof();
        for ( size_type i = 0; i < shiftScalar; ++i )
        {
            const size_type node = meshFEMFracture1.first_convex_of_basic_dof(i);

            nodo1 = meshFEMFracture2.point_of_basic_dof(node, 0);
            nodo2 = meshFEMFracture2.point_of_basic_dof(node, 1);

            const base_node midPoint = 0.5 * ( nodo1 + nodo2 );

            const scalar_type midPointSign = fractures[(f+1)%2]->getLevelSet()->getData()->levelSetFunction ( midPoint );

            for ( size_type kk = 0; kk < numberPoints; ++kk )
            {

                const base_node xx = 0.5 * (nodo1 + nodo2) + 0.5 * (nodo2 - nodo1)
                        * gauss_points [ kk ];

                const scalar_type gaussSign = fractures[(f+1)%2]->getLevelSet()->getData()->levelSetFunction ( xx );

                const base_node x0 = xx
                        - mesh->getMesh().points_of_convex(el) [ 0 ];

                const base_node x1 = xx
                        - mesh->getMesh().points_of_convex(el) [ 1 ];

                const base_node x2 = xx
                        - mesh->getMesh().points_of_convex(el) [ 2 ];

                const scalar_type A = 0.5 * (gmm::abs(x0 [ 0 ] * x1 [ 1 ]
                        - x0 [ 1 ] * x1 [ 0 ]) + gmm::abs(x0 [ 0 ] * x2 [ 1 ]
                        - x0 [ 1 ] * x2 [ 0 ]) + gmm::abs(x2 [ 0 ] * x1 [ 1 ]
                        - x2 [ 1 ] * x1 [ 0 ]));

                const scalar_type area = mesh->getMesh().convex_area_estimate(el);
                if ( gmm::abs(A - area) / area <= toll )
                {
                    if ( midPointSign * gaussSign > 0 )
                    {
                        (*(interpolationMatrices[f]))(j, i) += gauss_weights [ kk ];
                    }
                    continue;
                }
            }

        }

        for ( size_type ii = 0; ii < fractures[f]->getNumExtendedPressure(); ++ii )
        {
            const size_type i = fractures[f]->getExtendedPressure()[ii];
            const size_type node = meshFEMFracture1.first_convex_of_basic_dof(i);

            nodo1 = meshFEMFracture2.point_of_basic_dof ( node, 0);
            nodo2 = meshFEMFracture2.point_of_basic_dof ( node, 1);

            const base_node midPoint = 0.5 * ( nodo1 + nodo2 );

            const scalar_type midPointSign = fractures[(f+1)%2]->getLevelSet()->getData()->levelSetFunction ( midPoint );

            for ( size_type kk = 0; kk < numberPoints; ++kk )
            {

                const base_node xx = 0.5 * (nodo1 + nodo2) + 0.5 * (nodo2 - nodo1)
                        * gauss_points [ kk ];

                const scalar_type gaussSign = fractures[(f+1)%2]->getLevelSet()->getData()->levelSetFunction ( xx );

                const base_node x0 = xx
                        - mesh->getMesh().points_of_convex(el) [ 0 ];

                const base_node x1 = xx
                        - mesh->getMesh().points_of_convex(el) [ 1 ];

                const base_node x2 = xx
                        - mesh->getMesh().points_of_convex(el) [ 2 ];

                const scalar_type A = 0.5 * (gmm::abs(x0 [ 0 ] * x1 [ 1 ]
                        - x0 [ 1 ] * x1 [ 0 ]) + gmm::abs(x0 [ 0 ] * x2 [ 1 ]
                        - x0 [ 1 ] * x2 [ 0 ]) + gmm::abs(x2 [ 0 ] * x1 [ 1 ]
                        - x2 [ 1 ] * x1 [ 0 ]));

                const scalar_type area = mesh->getMesh().convex_area_estimate(el);
                if ( gmm::abs(A - area) / area <= toll )
                {
                    if ( midPointSign * gaussSign < 0 )
                    {
                        (*(interpolationMatrices[f]))(j, ii + fractures[f]->getMeshFEMPressure().nb_dof() ) += gauss_weights [ kk ];
                    }
                    continue;
                }
            }
        }

    }

    for ( size_type f = 0; f < fractures.size(); ++f )
    {
        scalar_type s=0;
        for ( size_type i = 0; i < fractures[f]->getMeshFEMPressure().nb_dof() + fractures[f]->getNumExtendedPressure(); ++i )
        {
            s += (*(interpolationMatrices[f]))(j, i);
        }

        if ( s > 0 )
        {
            for ( size_type i = 0; i < fractures[f]->getMeshFEMPressure().nb_dof() + fractures[f]->getNumExtendedPressure(); ++i )
            {
                (*(interpolationMatrices[f]))(j, i) = (*(interpolationMatrices[f]))(j, i) / s;
            }
        }
    }

} // interpolationMatrix


void massL2 ( sparseMatrixPtr_Type& M,
              const MeshHandlerPtr_Type& mesh,
              const scalarVector_Type& inv_eta,
              const size_type& uncutRegionFlag )
{

    const size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();

    generic_assembly assem;

    assem.set("w=data(#2);"
        "a=comp(vBase(#1).vBase(#1).Base(#2));"
        "M(#1,#1)+=a(:,i,:,i,j).w(j);");

    assem.push_mi(mesh->getIntegrationMethodScalar());

    assem.push_mf(mesh->getMeshFEMScalar());

    assem.push_mf(mesh->getMeshFEMScalar());

    assem.push_data(inv_eta);

    assem.push_mat(*M);

    assem.assembly(mesh->getRegion(uncutRegionFlag));

} // darcy_mass_p

// per il precondizionatore - e per la concentrazione, matrice di massa sulle p, per la frattura
void massL2 ( sparseMatrixPtr_Type& M,
              const FractureHandlerPtr_Type& fracture,
              const scalarVector_Type& dens )
{

    generic_assembly assem;

    if ( fracture->getMeshFEMPressure().get_qdim() == 1 )
    {
        assem.set("w=data(#2);"
            "a=comp(Base(#1).Base(#1).Base(#2));"
            "M(#1,#1)+=a(:, :, i).w(i);");
    }
    else
    {
        assem.set("w=data(#2);"
            "a=comp(vBase(#1).vBase(#1).Base(#2));"
            "M(#1,#1)+=a(:,i,:,i,j).w(j);");
    }

    assem.push_mi(fracture->getIntegrationMethodPressure());

    assem.push_mf(fracture->getMeshFEMPressure());
    assem.push_mf(fracture->getMeshFEMPressure());

    assem.push_data(dens);

    assem.push_mat(*M);

    assem.assembly(-1);

} // darcy_mass_pF

void massL2 ( sparseMatrixPtr_Type& M,
              const MeshHandlerPtr_Type& mesh,
              const FractureHandlerPtr_Type& fracture,
              const scalarVector_Type& inv_eta,
              const size_type& cutRegionFlag )
{
    const scalar_type fractureID = fracture->getId();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const sizeVector_Type& extendedDOF = mesh->getExtendedDOFScalar(fractureID);
    const size_type shiftExtended = extendedDOF.size();
    size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();

    sparseMatrix_Type MIn, MOut;
    gmm::resize(MIn, shiftPressure, shiftPressure);
    gmm::resize(MOut, shiftPressure, shiftPressure);
    gmm::clear(MIn);
    gmm::clear(MOut);

    shiftPressure += mesh->getCountExtendedDOFScalar(fractureID - 1);
    ;

    generic_assembly assemIn, assemOut;

    assemIn.set("w=data(#2);"
        "a=comp(vBase(#1).vBase(#1).Base(#2));"
        "M(#1,#1)+=a(:,i,:,i,j).w(j);");

    assemOut.set("w=data(#2);"
        "a=comp(vBase(#1).vBase(#1).Base(#2));"
        "M(#1,#1)+=a(:,i,:,i,j).w(j);");

    assemIn.push_mi(levelSet->getIntegrationMethodInside());
    assemOut.push_mi(levelSet->getIntegrationMethodOutside());

    assemIn.push_mf(mesh->getMeshFEMScalar());
    assemOut.push_mf(mesh->getMeshFEMScalar());

    assemIn.push_mf(mesh->getMeshFEMScalar());
    assemOut.push_mf(mesh->getMeshFEMScalar());

    assemIn.push_data(inv_eta);
    assemOut.push_data(inv_eta);

    assemIn.push_mat(MIn);
    assemOut.push_mat(MOut);

    assemIn.assembly(mesh->getRegion(cutRegionFlag));
    assemOut.assembly(mesh->getRegion(cutRegionFlag));

    //qui devo aggiungere i gdl raddoppiati
    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = extendedDOF [ i ];
        for ( size_type j = 0; j < shiftExtended; ++j )
        {
            size_type jj = extendedDOF [ j ];
            if ( (levelSet->getBaricenterValue(ii) <= 0)
                    && (levelSet->getBaricenterValue(jj) <= 0) )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn(ii, jj);
                (*M)(i + shiftPressure, j + shiftPressure) += MOut(ii, jj);
            }
            else if ( (levelSet->getBaricenterValue(ii) >= 0)
                    && (levelSet->getBaricenterValue(jj) >= 0) )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + shiftPressure, j + shiftPressure) += MIn(ii, jj);
            }
            else if ( (levelSet->getBaricenterValue(ii) <= 0)
                    && (levelSet->getBaricenterValue(jj) >= 0) )
            {
                // i is In, j is Out
                (*M)(ii, j + shiftPressure) += MIn(ii, jj);
                (*M)(i + shiftPressure, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + shiftPressure, jj) += MIn(ii, jj);
                (*M)(ii, j + shiftPressure) += MOut(ii, jj);
            }
        }
    }
} // darcy_mass_p

void massHdiv ( sparseMatrixPtr_Type& M,
                const MeshHandlerPtr_Type& mesh,
                const scalarVectorContainer_Type& mediumEtaInterpolated,
                const size_type& uncutRegionFlag )
{

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize(M_, shiftVelocity, shiftVelocity);
    gmm::clear(M_);

    //questo è il termine tau tau (matrice di "massa" per le velocità)

    generic_assembly assem;

//    assem.set("w=data(#2);"
//        "a=comp(vBase(#1).vBase(#1).Base(#2));"
//        "M(#1,#1)+=a(:,i,:,i,j).w(j);");

    assem.set ( "xx=data(#2);""yy=data$2(#2);""xy=data$3(#2);"
                "a=comp(vBase(#1).vBase(#1).Base(#2));"
                "M(#1,#1)+=a(:,1,:,1,j).xx(j)+a(:,2,:,2,j).yy(j)+a(:,1,:,2,j).xy(j)+a(:,2,:,1,j).xy(j);");

    // Assign the mesh integration method
    assem.push_mi(mesh->getIntegrationMethodVector());

    // Assign the mesh finite element space
    assem.push_mf(mesh->getMeshFEMVector());

    // Assign the mesh finite element space for the coefficients
    assem.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the coefficients
    assem.push_data(mediumEtaInterpolated[0]);

    assem.push_data(mediumEtaInterpolated[2]);

    assem.push_data(mediumEtaInterpolated[1]);

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly(mesh->getRegion(uncutRegionFlag));

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    }

    cout << "DARCY :: operator a(volume)      [OK]" << endl;

}

void massHdiv ( sparseMatrixPtr_Type& M,
                const MeshHandlerPtr_Type& mesh,
                const FractureHandlerPtr_Type& fracture,
                const scalarVector_Type& fractureEtaNormalOnMedium,
                const scalarVectorContainer_Type& mediumEtaInterpolated,
                const size_type& cutRegionFlag )
{

    const scalar_type fractureID = fracture->getId();
    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const sizeVector_Type extendedDOF = mesh->getExtendedDOFVector(fractureID);
    const size_type shiftExtended = extendedDOF.size();

    sparseMatrix_Type MIn, MOut, Gamma;
    gmm::resize(MOut, shiftVelocity, shiftVelocity);
    gmm::clear(MOut);
    gmm::resize(MIn, shiftVelocity, shiftVelocity);
    gmm::clear(MIn);
    gmm::resize(Gamma, shiftVelocity, shiftVelocity);
    gmm::clear(Gamma);

    shiftVelocity += mesh->getCountExtendedDOFVector(fractureID - 1);

    //questo è il termine tau tau (matrice di "massa" per le velocità)

    generic_assembly assemIn, assemOut, assemGam;

//    assemIn.set("w=data(#2);"
//        "a=comp(vBase(#1).vBase(#1).Base(#2));"
//        "M(#1,#1)+=a(:,i,:,i,j).w(j);");

//    assemOut.set("w=data(#2);"
//        "a=comp(vBase(#1).vBase(#1).Base(#2));"
//        "M(#1,#1)+=a(:,i,:,i,j).w(j);");

    assemIn.set ( "xx=data(#2);""yy=data$2(#2);""xy=data$3(#2);"
                  "a=comp(vBase(#1).vBase(#1).Base(#2));"
                  "M(#1,#1)+=a(:,1,:,1,j).xx(j)+a(:,2,:,2,j).yy(j)+a(:,1,:,2,j).xy(j)+a(:,2,:,1,j).xy(j);" );

    assemOut.set ( "xx=data(#2);""yy=data$2(#2);""xy=data$3(#2);"
                   "a=comp(vBase(#1).vBase(#1).Base(#2));"
                   "M(#1,#1)+=a(:,1,:,1,j).xx(j)+a(:,2,:,2,j).yy(j)+a(:,1,:,2,j).xy(j)+a(:,2,:,1,j).xy(j);" );

    // Assign the mesh integration method
    assemIn.push_mi(levelSet->getIntegrationMethodInside());
    assemOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the mesh finite element space
    assemIn.push_mf(mesh->getMeshFEMVector());
    assemOut.push_mf(mesh->getMeshFEMVector());

    // Assign the mesh finite element space for the coefficients
    assemIn.push_mf(mesh->getMeshFEMCoefficients());
    assemOut.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the coefficients
    assemIn.push_data(mediumEtaInterpolated[0]);
    assemIn.push_data(mediumEtaInterpolated[2]);
    assemIn.push_data(mediumEtaInterpolated[1]);

    assemOut.push_data(mediumEtaInterpolated[0]);
    assemOut.push_data(mediumEtaInterpolated[2]);
    assemOut.push_data(mediumEtaInterpolated[1]);

    // Set the matrices to save the evaluations
    assemIn.push_mat(MIn);
    assemOut.push_mat(MOut);

    // Computes the matrices
    assemIn.assembly(mesh->getRegion(cutRegionFlag));
    assemOut.assembly(mesh->getRegion(cutRegionFlag));

    cout << "DARCY :: operator a(volume)      [OK]" << endl;

    // Add the extended degrees of freedom
    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = extendedDOF [ i ];

        for ( size_type j = 0; j < shiftExtended; ++j )
        {
            size_type jj = extendedDOF [ j ];

            if ( (levelSet->getDOFValue(ii) < 0) && (levelSet->getDOFValue(jj)
                    < 0) )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) += MOut(ii, jj);
            }
            else
            {
                if ( (levelSet->getDOFValue(ii) >= 0)
                        && (levelSet->getDOFValue(jj) >= 0) )
                {
                    // i and j are both Out
                    (*M)(ii, jj) += MOut(ii, jj);
                    (*M)(i + shiftVelocity, j + shiftVelocity) += MIn(ii, jj);
                }
                else
                {
                    if ( (levelSet->getDOFValue(ii) < 0)
                            && (levelSet->getDOFValue(jj) >= 0) )
                    {
                        // i is In, j is Out
                        (*M)(ii, j + shiftVelocity) += MIn(ii, jj);
                        (*M)(i + shiftVelocity, jj) += MOut(ii, jj);
                    }
                    else
                    {
                        // i is Out, j is In
                        (*M)(i + shiftVelocity, jj) += MIn(ii, jj);
                        (*M)(ii, j + shiftVelocity) += MOut(ii, jj);
                    }
                }
            }
        }


    }
    // Computes the integral {u dot n}{v dot n}

    assemGam.set("w=data(#2);"
        "a=comp(vBase(#1).NonLin(#3).vBase(#1).NonLin(#3).Base(#2));"
        "M(#1,#1)+=a(:,j,j,:,i,i,k).w(k);");

    level_set_unit_normal nterm(levelSet->getLevelSet().get_mesh_fem(),
            levelSet->getLevelSet().values());

    // Assign the mesh integration method
    assemGam.push_mi(levelSet->getIntegrationMethod());

    // Assign the mesh finite element space
    assemGam.push_mf(mesh->getMeshFEMVector());
    assemGam.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the non linear term
    assemGam.push_nonlinear_term(&nterm);

    // Assign the mesh finite element space for the coefficients
    assemGam.push_mf(levelSet->getLevelSet().get_mesh_fem());

    // Assign the coefficients
    assemGam.push_data(fractureEtaNormalOnMedium);

    // Set the matrices to save the evaluations
    assemGam.push_mat(Gamma);

    // Computes the matrices
    assemGam.assembly(mesh->getRegion(cutRegionFlag));

    // Add the extended degrees of freedom
    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = mesh->getExtendedDOFVector(fractureID, i);
        for ( size_type j = 0; j < shiftExtended; ++j )
        {
            size_type jj = mesh->getExtendedDOFVector(fractureID, j);
            (*M)(ii, jj) += 0.25 * Gamma(ii, jj);
            (*M)(i + shiftVelocity, j + shiftVelocity) += 0.25 * Gamma(ii, jj);
            (*M)(i + shiftVelocity, jj) += 0.25 * Gamma(ii, jj);
            (*M)(ii, j + shiftVelocity) += 0.25 * Gamma(ii, jj);
        }
    }

    // Computes the integral [u dot n][v dot n].
    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = mesh->getExtendedDOFVector(fractureID, i);
        for ( size_type j = 0; j < shiftExtended; ++j )
        {
            size_type jj = mesh->getExtendedDOFVector(fractureID, j);

            // If v (ii) is In, vOut = 0, the jump is positive
            if ( levelSet->getDOFValue(ii) * levelSet->getDOFValue(jj) >= 0 )
            {
                (*M)(ii, jj) += 0.5 * fracture->getData().getCsi0() * Gamma(ii,
                        jj);
                (*M)(i + shiftVelocity, jj) -= 0.5
                        * fracture->getData().getCsi0() * Gamma(ii, jj);
                (*M)(ii, j + shiftVelocity) -= 0.5
                        * fracture->getData().getCsi0() * Gamma(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) += 0.5
                        * fracture->getData().getCsi0() * Gamma(ii, jj);

            }

            // If v (ii) is Out, vIn = 0, the jump is negative

            else
            {
                (*M)(ii, jj) -= 0.5 * fracture->getData().getCsi0() * Gamma(ii,
                        jj);
                (*M)(i + shiftVelocity, jj) += 0.5
                        * fracture->getData().getCsi0() * Gamma(ii, jj);
                (*M)(ii, j + shiftVelocity) += 0.5
                        * fracture->getData().getCsi0() * Gamma(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) -= 0.5
                        * fracture->getData().getCsi0() * Gamma(ii, jj);

            }
        }
    }

}

void massHdiv ( sparseMatrixPtr_Type& M, const IntersectData_Type& intersect, const MeshHandlerPtr_Type& mesh,
                const scalarVectorContainer_Type& mediumEtaInterpolated,
                const std::vector < scalarVector_Type* >& fracturesEta )
{
    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    sparseMatrix_Type Matrix;
    gmm::resize(Matrix, shiftVelocity, shiftVelocity);

    getfem::mesh_fem::ind_dof_ct indexDOF = mesh->getMeshFEMVector().ind_basic_dof_of_element ( intersect.getElementID() );
    const size_type indexDOFsize = indexDOF.size();
    const size_type regionIDsize = intersect.getRegionActive().size();

    for ( size_type regionID = 0; regionID < regionIDsize; ++regionID )
    {
        gmm::clear ( Matrix );
        generic_assembly assembly;
        //assembly.set("w=data(#2);" "a=comp(vBase(#1).vBase(#1).Base(#2));" "M(#1,#1)+=a(:,i,:,i,j).w(j);");
        assembly.set ( "xx=data(#2);""yy=data$2(#2);""xy=data$3(#2);"
                       "a=comp(vBase(#1).vBase(#1).Base(#2));"
                       "M(#1,#1)+=a(:,1,:,1,j).xx(j)+a(:,2,:,2,j).yy(j)+a(:,1,:,2,j).xy(j)+a(:,2,:,1,j).xy(j);" );

        // construct the integration method based on the meshLevelSet
        const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_TRIANGLE(6)" );
        getfem::mesh_im_level_set meshImLevelSet ( mesh->getMeshLevelSet(), getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );
        meshImLevelSet.set_integration_method ( mesh->getMesh().convex_index(), intTypeIM );
        meshImLevelSet.set_simplex_im ( intTypeIM );

        getfem::mesh_region meshElement;
        meshElement.add ( intersect.getElementID() );

        sizeVector_Type levelSetID ( intersect.getNumFractures(), 0 );
        for ( size_type f = 0; f < levelSetID.size(); ++f )
        {
            levelSetID [f] = intersect.getFracture(f)->getId();
        }

        assembly.push_mi ( meshImLevelSet );
        assembly.push_mf ( mesh->getMeshFEMVector() );

        // Assign the mesh finite element space for the coefficients
        assembly.push_mf(mesh->getMeshFEMCoefficients());

        // Assign the coefficients
        assembly.push_data(mediumEtaInterpolated[0]);
        assembly.push_data(mediumEtaInterpolated[2]);
        assembly.push_data(mediumEtaInterpolated[1]);

        assembly.push_mat ( Matrix );

        // Set the operation bewteen the level sets
        const std::string activeRegion = intersect.getRegionActive()[regionID];

        const std::string operation = getOperation ( activeRegion, levelSetID );

        meshImLevelSet.set_level_set_boolean_operations ( operation );

        meshImLevelSet.adapt();

        assembly.assembly ( meshElement );

        // save the solution in M
        // loop on rows
        for ( size_type i = 0; i < indexDOFsize; ++i )
        {
            const size_type shiftI = intersect.getDOFVelocity ( i, regionID );

            // loop on column
            for ( size_type j = 0; j < indexDOFsize; ++j )
            {
                const size_type shiftJ = intersect.getDOFVelocity ( j, regionID );
                (*M)( shiftI, shiftJ ) += Matrix ( indexDOF[i], indexDOF[j] );
            }
        }
    }

    const size_type numFractures = intersect.getFractures().size();

    // Computes the integral {u dot n}{v dot n} and [u dot n] [v dot n]
    for ( size_type f = 0; f < numFractures; ++f )
    {
        const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_TRIANGLE(6)" );
        LevelSetHandlerPtr_Type& levelSetFracture = intersect.getFracture(f)->getLevelSet();
        getfem::level_set levelSet ( mesh->getMesh(), dim_type(1), true );

        levelSet.values(0) = levelSetFracture->getLevelSet().values(0);
        scalarVector_Type levelSetvalues1 = intersect.getFracture ( ((int)f + 1)%2 )->getLevelSet()->getLevelSet().values(0);

        sizeVector_Type facingRegion = intersect.getFacingRegion ( f );

        for ( size_type segment = 0; segment < 2; ++segment )
        {
            for ( size_type ii = 0; ii < levelSetvalues1.size(); ++ii )
            {
                levelSetvalues1 [ ii ] *= std::pow( -1, segment );
            }

             levelSet.values(1) = levelSetvalues1;

            getfem::mesh_level_set meshLevelSet ( mesh->getMesh() );
            meshLevelSet.add_level_set ( levelSet );
            meshLevelSet.adapt();

            getfem::mesh_im_level_set imLevelSet ( meshLevelSet, getfem::mesh_im_level_set::INTEGRATE_BOUNDARY );
            imLevelSet.set_integration_method ( intersect.getElementID(), intTypeIM );
            imLevelSet.set_simplex_im ( intTypeIM );

            imLevelSet.adapt();

            generic_assembly assembly;

            assembly.set ( "w=data(#2);"
                           "a=comp(vBase(#1).NonLin(#3).vBase(#1).NonLin(#3).Base(#2));"
                           "M(#1,#1)+=a(:,j,j,:,i,i,k).w(k);");

            level_set_unit_normal nterm ( levelSetFracture->getLevelSet().get_mesh_fem(),
                                          levelSetFracture->getLevelSet().values() );

            getfem::mesh_region meshElement;
            meshElement.add ( intersect.getElementID() );

            // Assign the mesh integration method
            assembly.push_mi ( imLevelSet );

            // Assign the mesh finite element space
            assembly.push_mf ( mesh->getMeshFEMVector() );
            assembly.push_mf ( mesh->getMeshFEMCoefficients() );

            // Assign the non linear term
            assembly.push_nonlinear_term ( &nterm );

            // Assign the mesh finite element space for the coefficients
            assembly.push_mf ( levelSetFracture->getLevelSet().get_mesh_fem() );

            // Assign the coefficients
            assembly.push_data ( *(fracturesEta[f]) );

            // Set the matrices to save the evaluations
            sparseMatrix_Type Gamma;
            const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
            gmm::resize(Gamma, shiftVelocity, shiftVelocity);
            gmm::clear(Gamma);
            assembly.push_mat ( Gamma );

            // Computes the matrices
            assembly.assembly ( meshElement );

            getfem::mesh_fem::ind_dof_ct indexDOFVector = mesh->getMeshFEMVector().ind_basic_dof_of_element ( intersect.getElementID() );

            // Add the extended degrees of freedom
            for ( size_type regionID1 = 0; regionID1 < 2; ++regionID1 )
            {
                const size_type index1= regionID1 + segment * ( facingRegion.size() - 2 );
                for ( size_type regionID2 = 0; regionID2 < 2; ++regionID2 )
                {
                    const size_type index2 = regionID2 + segment * ( facingRegion.size() - 2 );
                    for ( size_type i = 0; i < indexDOFVector.size(); ++i )
                    {
                        const size_type shiftI = intersect.getDOFVelocity ( i, facingRegion[index1] );
                        for ( size_type j = 0; j < indexDOFVector.size(); ++j )
                        {
                            const size_type shiftJ = intersect.getDOFVelocity ( j, facingRegion[index2] );
                            (*M)( shiftI, shiftJ) += 0.25 * Gamma ( indexDOFVector[i], indexDOFVector[j] );
                            if ( regionID1 == regionID2 )
                            {
                                (*M)( shiftI, shiftJ) += 0.5 * intersect.getFracture(f)->getData().getCsi0() *
                                                        Gamma ( indexDOFVector[i], indexDOFVector[j] );
                            }
                            else
                            {
                                (*M)( shiftI, shiftJ) -= 0.5 * intersect.getFracture(f)->getData().getCsi0() *
                                                        Gamma ( indexDOFVector[i], indexDOFVector[j] );
                            }
                        }
                    }
                }
            }

        }//fine for sui segment
    }

} // massHdiv

void darcy_E ( sparseMatrixPtr_Type& E,
               const MeshHandlerPtr_Type& mesh,
               const FractureHandlerPtr_Type& fracture,
               const size_type& cutRegionFlag )
{
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const scalar_type fractureID = fracture->getId();
    const sizeVector_Type& extendedDOFScalar = mesh->getExtendedDOFScalar(
            fractureID);
    const size_type shiftExtendedScalar = extendedDOFScalar.size();
    const sizeVector_Type& extendedDOFVector = mesh->getExtendedDOFVector(
            fractureID);
    const size_type shiftExtendedVector = extendedDOFVector.size();

    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();

    sparseMatrix_Type F;
    gmm::resize(F, shiftPressure, shiftVelocity);
    gmm::clear(F);

    shiftPressure += mesh->getCountExtendedDOFScalar(fractureID - 1);
    shiftVelocity += mesh->getCountExtendedDOFVector(fractureID - 1);

    generic_assembly assemF;

    // Computes the integral {u dot n}q. Term on Gamma for the join of the solution in the fracure and outside
    // sarebbe {u dot n}pf e il suo trasposto serve per [u dot n]phi (termine sorgente nella frattura)

    assemF.set("a=comp(Base(#2).vBase(#1).NonLin(#3));"
        "M(#2,#1)+=a(:,:,j,j);");

    // Assign the mesh integration method
    assemF.push_mi(levelSet->getIntegrationMethod());

    // Assign the mesh finite element space
    assemF.push_mf(mesh->getMeshFEMVector());
    assemF.push_mf(mesh->getMeshFEMScalar());
    assemF.push_mf(levelSet->getLevelSet().get_mesh_fem());

    // Assign the non linear term
    level_set_unit_normal nterm(levelSet->getLevelSet().get_mesh_fem(),
            levelSet->getLevelSet().values());
    assemF.push_nonlinear_term(&nterm);

    // Set the matrices to save the evaluations
    assemF.push_mat(F);

    // Computes the matrices
    assemF.assembly(mesh->getRegion(cutRegionFlag));

    // Update the extended degrees of freedom
    for ( size_type i = 0; i < shiftExtendedVector; ++i )
    {
        size_type ii = extendedDOFVector [ i ];
        for ( size_type j = 0; j < shiftExtendedScalar; ++j )
        {
            size_type jj = extendedDOFScalar [ j ];

            if ( levelSet->getDOFValue(ii) < 0 )
            {
                (*E)(ii, jj) += F(jj, ii);
                (*E)(i + shiftVelocity, jj) -= F(jj, ii);
            }
            else
            {
                (*E)(ii, jj) -= F(jj, ii);
                (*E)(i + shiftVelocity, jj) += F(jj, ii);
            }

        }
    }

} // darcy_E

void darcy_E ( sparseMatrixPtr_Type& E,
               const FractureHandlerPtr_Type& fracture,
               const IntersectData_Type& intersect,
               const MeshHandlerPtr_Type& mesh,
               const scalar_type& localFractureID )
{

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_TRIANGLE(6)" );
    LevelSetHandlerPtr_Type& levelSetFracture = fracture->getLevelSet();
    getfem::level_set levelSet ( mesh->getMesh(), dim_type(1), true );

    levelSet.values(0) = levelSetFracture->getLevelSet().values(0);
    scalarVector_Type levelSetvalues1 = intersect.getFracture ( ((int)localFractureID + 1)%2 )->getLevelSet()->getLevelSet().values(0);

    sizeVector_Type facingRegion = intersect.getFacingRegion ( localFractureID );

    for ( size_type segment = 0; segment < intersect.getFractures().size(); ++segment )
    {
        for ( size_type ii = 0; ii < levelSetvalues1.size(); ++ii )
        {
            levelSetvalues1 [ ii ] *= std::pow( -1, segment );
        }

        levelSet.values(1) = levelSetvalues1;

        getfem::mesh_level_set meshLevelSet ( mesh->getMesh() );
        meshLevelSet.add_level_set ( levelSet );
        meshLevelSet.adapt();

        getfem::mesh_im_level_set imLevelSet ( meshLevelSet, getfem::mesh_im_level_set::INTEGRATE_BOUNDARY );
        imLevelSet.set_integration_method ( intersect.getElementID(), intTypeIM );
        imLevelSet.set_simplex_im ( intTypeIM );

        imLevelSet.adapt();

        generic_assembly assemF;
        assemF.set("a=comp(Base(#2).vBase(#1).NonLin(#3));"
                   "M(#2,#1)+=a(:,:,j,j);");
        assemF.push_mi ( imLevelSet );

        getfem::mesh_region meshElement;
        meshElement.add ( intersect.getElementID() );

        // Assign the mesh finite element space
        assemF.push_mf ( mesh->getMeshFEMVector() );
        assemF.push_mf ( mesh->getMeshFEMScalar() );
        assemF.push_mf ( levelSetFracture->getLevelSet().get_mesh_fem() );

        // Assign the non linear term
        level_set_unit_normal nterm( levelSetFracture->getLevelSet().get_mesh_fem(),
                                     levelSetFracture->getLevelSet().values() );
        assemF.push_nonlinear_term ( &nterm );

        // Set the matrices to save the evaluations
        size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
        size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();

        sparseMatrix_Type F;
        gmm::resize(F, shiftPressure, shiftVelocity);
        gmm::clear(F);

        assemF.push_mat(F);

        assemF.assembly ( meshElement );

        getfem::mesh_fem::ind_dof_ct indexDOFVector = mesh->getMeshFEMVector().ind_basic_dof_of_element ( intersect.getElementID() );
        getfem::mesh_fem::ind_dof_ct indexDOFScalar = mesh->getMeshFEMScalar().ind_basic_dof_of_element ( intersect.getElementID() );

        // loop on rows
        scalarVector_Type sign (2);
        sign[0] = 1; sign[1] = -1;

        for ( size_type regionID = 0; regionID < 2; ++regionID )
        {
            for ( size_type i = 0; i < indexDOFVector.size(); ++i )
            {
                const size_type index = regionID + segment * ( facingRegion.size() - 2 );
                const size_type shiftI = intersect.getDOFVelocity ( i, facingRegion[index] );
                (*E)( shiftI, indexDOFScalar [ 0 ] ) += -sign[(regionID+1)%2] * F ( indexDOFScalar [ 0 ], indexDOFVector [ i ] );
            }
        }
    }

} // darcy_E

void nitsche ( sparseMatrixPtr_Type& M,
               const MeshHandlerPtr_Type& mesh,
               const scalarVectorContainer_Type& mediumInvKInterpolated,
               const scalar_type& penaltyParameterVelocity,
               const sizeVector_Type& extBoundaryUncut,
               const size_type& uncutRegionFlag )
{
    // Boundary conditions, usign the Nitche penalization

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    const size_type shiftCoefficients = mesh->getMeshFEMCoefficients().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize(M_, shiftVelocity, shiftVelocity);
    gmm::clear(M_);

    scalarVector_Type etaGammaUinvh(shiftCoefficients);

    for ( size_type i = 0; i < shiftCoefficients; ++i )
    {
        const scalar_type firstRow = mediumInvKInterpolated[0] [ i ] + mediumInvKInterpolated[1] [ i ];
        const scalar_type secondRow = mediumInvKInterpolated[2] [ i ] + mediumInvKInterpolated[1] [ i ];
        const scalar_type norminf = std::max ( firstRow, secondRow );
        etaGammaUinvh [ i ] = norminf * penaltyParameterVelocity * mesh->getInverseMeshSizeDOF(i);
    }

    getfem::generic_assembly assem_surf;

    assem_surf.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(mesh->getIntegrationMethodVector());

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(mesh->getMeshFEMVector());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the coefficients
    assem_surf.push_data(etaGammaUinvh);

    // Set the matrices to save the evaluations
    assem_surf.push_mat(M_);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < extBoundaryUncut.size(); bndID++ )
    {
        assem_surf.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(
                        extBoundaryUncut [ bndID ]));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) += M_(i, j);
        }
    }

    cout << "DARCY :: operator a(surface)     [OK]" << endl;

} // nitsche

// Boundary conditions, usign the Nitche penalization
void nitsche ( sparseMatrixPtr_Type& M,
               const MeshHandlerPtr_Type& mesh,
               const FractureHandlerPtr_Type& fracture,
               const scalarVectorContainer_Type& mediumInvKInterpolated,
               const scalar_type& penaltyParameterVelocity,
               const sizeVector_Type& ExtBoundary_cut )
{
    const scalar_type fractureID = fracture->getId();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const sizeVector_Type& extendedDOF = mesh->getExtendedDOFVector(fractureID);
    const size_type shiftExtended = extendedDOF.size();
    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    const size_type shiftCoefficients = mesh->getMeshFEMCoefficients().nb_dof();

    sparseMatrix_Type MIn, MOut;
    gmm::resize(MOut, shiftVelocity, shiftVelocity);
    gmm::clear(MOut);
    gmm::resize(MIn, shiftVelocity, shiftVelocity);
    gmm::clear(MIn);

    shiftVelocity += mesh->getCountExtendedDOFVector(fractureID - 1);

    scalarVector_Type etaGammaUinvh(shiftCoefficients);

    for ( size_type i = 0; i < shiftCoefficients; ++i )
    {
        const scalar_type firstRow = mediumInvKInterpolated[0] [ i ] + mediumInvKInterpolated[1] [ i ];
        const scalar_type secondRow = mediumInvKInterpolated[2] [ i ] + mediumInvKInterpolated[1] [ i ];
        const scalar_type norminf = std::max ( firstRow, secondRow );
        etaGammaUinvh [ i ] = norminf * penaltyParameterVelocity * mesh->getInverseMeshSizeDOF(i);
    }

    getfem::generic_assembly assem_surfIn, assem_surfOut;

    assem_surfIn.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    assem_surfOut.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assign the mesh integration method
    assem_surfIn.push_mi(levelSet->getIntegrationMethodInside());
    assem_surfOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the mesh finite element space
    assem_surfIn.push_mf(mesh->getMeshFEMVector());
    assem_surfOut.push_mf(mesh->getMeshFEMVector());

    // Assign the mesh finite element space for the coefficients
    assem_surfIn.push_mf(mesh->getMeshFEMCoefficients());
    assem_surfOut.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the coefficients
    assem_surfIn.push_data(etaGammaUinvh);
    assem_surfOut.push_data(etaGammaUinvh);

    // Set the matrices to save the evaluations
    assem_surfIn.push_mat(MIn);
    assem_surfOut.push_mat(MOut);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < ExtBoundary_cut.size(); bndID++ )
    {
        assem_surfIn.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(
                        ExtBoundary_cut [ bndID ]));
        assem_surfOut.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(
                        ExtBoundary_cut [ bndID ]));
    }
    // Update the extended degrees of freedom
    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = extendedDOF [ i ];
        for ( size_type j = 0; j < shiftExtended; ++j )
        {
            size_type jj = extendedDOF [ j ];
            if ( (levelSet->getDOFValue(ii) < 0) && (levelSet->getDOFValue(jj)
                    < 0) )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) += MOut(ii, jj);
            }

            else if ( (levelSet->getDOFValue(ii) >= 0)
                    && (levelSet->getDOFValue(jj) >= 0) )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) += MIn(ii, jj);
            }
            else if ( (levelSet->getDOFValue(ii) < 0)
                    && (levelSet->getDOFValue(jj) >= 0) )
            {
                // i is In, j is Out
                (*M)(ii, j + shiftVelocity) += MIn(ii, jj);
                (*M)(i + shiftVelocity, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + shiftVelocity, jj) += MIn(ii, jj);
                (*M)(ii, j + shiftVelocity) += MOut(ii, jj);
            }
        }
    }

    cout << "DARCY :: operator a(surface)     [OK]" << endl;

} // nitsche

//matrice tau tau per la frattura
void darcy_A11F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const scalar_type& gammaU,
                  const scalarVector_Type& invKTangentialInterpolated,
                  const sizeVector_Type& ExtBoundary,
                  const size_type& uncutRegionFlag )
{
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    const size_type shiftData = fracture->getMeshFEMPressure().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize(M_, shiftVelocity, shiftVelocity);
    gmm::clear(M_);

    scalarVector_Type etaGammaUinvh(shiftData);

    for ( size_type i = 0; i < shiftData; ++i )
    {
        etaGammaUinvh [ i ] = invKTangentialInterpolated [ i ] * gammaU
                * fracture->getInverseMeshSize(i);
    }

    // Volume integration
    const size_type shiftMapFactor =
            fracture->getMagnificationMapFactor1().size();
    scalarVector_Type invF(shiftMapFactor, 0.);

    generic_assembly assem, assemGam;
    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        assem.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
            "M(#1,#1)+=a(:, :, i,k).w(i).q(k);");
    }
    else
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1 / (fracture->getMagnificationMapFactor1(i)
                    * fracture->getMagnificationMapFactor2(i));
        }
        assem.set("w=data(#2);" "q=data$2(#2);"
            "a=comp(vBase(#1).vBase(#1).Base(#2).Base(#2));"
            "M(#1,#1)+=a(:,i,:,i,j,k).w(j).q(k);");
    }

    // Assign the M_mediumMesh integration method
    assem.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem.push_mf(fracture->getMeshFEMPressure());
    assem.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem.push_data(invKTangentialInterpolated);
    assem.push_data(invF);

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly ( uncutRegionFlag );

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operator a(volume)      [OK]" << endl;

    // Boundary integration for the fracture
    gmm::clear(M_);

    getfem::generic_assembly assem_surf, assem_normx, assem_normy;

    assem_surf.set("gamma=data$1(#2);"
        "t=comp(vBase(#1).Normal().vBase(#1).Normal().Base(#2));"
        "M$1(#1,#1)+=(t(:,i, i, :,j, j, k).gamma(k));");

    // Assign the M_mediumMesh integration method
    assem_surf.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem_surf.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_surf.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem_surf.push_data(etaGammaUinvh);

    // Set the matrices to save the evaluations
    assem_surf.push_mat(M_);

    // Assemble in each sub region
    for ( size_type bndID = 0; bndID < ExtBoundary.size(); bndID++ )
    {
        assem_surf.assembly(
                fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region(
                        ExtBoundary [ bndID ]));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftVelocity; ++j )
        {
            (*M)(i, j) += M_(i, j);
        }
    }
    cout << "DARCY :: operator a(surface)     [OK]" << endl;

} // darcy_A11F

//matrice tau tau per la frattura intersecata
void darcy_A11F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const scalarVector_Type& invKTangentialInterpolated,
                  const FractureHandlerPtr_Type& otherFracture,
                  const size_type& cutRegionFlag )
{

    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    sparseMatrix_Type MIn, MOut, Gamma;
    gmm::resize(MOut, shiftVelocity, shiftVelocity);
    gmm::clear(MOut);
    gmm::resize(MIn, shiftVelocity, shiftVelocity);
    gmm::clear(MIn);
    gmm::resize(Gamma, shiftVelocity, shiftVelocity);
    gmm::clear(Gamma);

    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    const size_type shiftMapFactor =
                    fracture->getMagnificationMapFactor1().size();
    scalarVector_Type invF(shiftMapFactor, 0.);

    for ( size_type i = 0; i < shiftMapFactor; ++i )
    {
        invF [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
    }

    generic_assembly assemIn, assemOut, assemGam;

    assemIn.set("w=data(#2);" "q=data$2(#2);"
                "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
                "M(#1,#1)+=a(:,:,j,k).w(j).q(k);");

    assemOut.set("w=data(#2);" "q=data$2(#2);"
                 "a=comp(Base(#1).Base(#1).Base(#2).Base(#2));"
                 "M(#1,#1)+=a(:,:,j,k).w(j).q(k);");

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevelSetOut ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                  getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );

    getfem::mesh_im_level_set meshImLevelSetIn ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                 getfem::mesh_im_level_set::INTEGRATE_INSIDE );

    meshImLevelSetOut.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );
    meshImLevelSetIn.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevelSetOut.set_simplex_im ( intTypeIM );
    meshImLevelSetIn.set_simplex_im ( intTypeIM );


    // Assign the mesh integration method
    assemIn.push_mi ( meshImLevelSetIn );
    assemOut.push_mi ( meshImLevelSetOut );

    // Assign the mesh finite element space
    assemIn.push_mf ( fracture->getMeshFEMVelocity() );
    assemOut.push_mf ( fracture->getMeshFEMVelocity() );

    // Assign the mesh finite element space for the coefficients
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Assign the coefficients
    assemIn.push_data ( invKTangentialInterpolated );
    assemIn.push_data ( invF );

    assemOut.push_data ( invKTangentialInterpolated );
    assemOut.push_data ( invF );

    assemOut.push_mat ( MOut );
    assemIn.push_mat ( MIn );

    assemOut.assembly ( cutRegionFlag );
    assemIn.assembly ( cutRegionFlag );

    // Update the extended degrees of freedom
    const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();
    const size_type extendedNumVelocity = fracture->getNumExtendedVelocity();

    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
        const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
        base_node pointMapped(0,0);
        pointMapped[0] = pointFlat[0];
        pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
        const scalar_type levelSetValue1 = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );

        for ( size_type j = 0; j < extendedNumVelocity; ++j )
        {
            const size_type jj = extendedVelocity [ j ];
            const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(jj);
            base_node pointMapped(0,0);
            pointMapped[0] = pointFlat[0];
            pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
            const scalar_type levelSetValue2 = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );

            if ( levelSetValue1 < 0 && levelSetValue2 < 0 )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn ( ii, jj );
                (*M)(i + shiftVelocity, j + shiftVelocity) += MOut ( ii, jj );
            }

            else if ( levelSetValue1 >= 0 && levelSetValue2 >= 0 )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + shiftVelocity, j + shiftVelocity) += MIn(ii, jj);
            }
            else if ( levelSetValue1 < 0 && levelSetValue2 >= 0 )
            {
                // i is In, j is Out
                (*M)(ii, j + shiftVelocity) += MIn(ii, jj);
                (*M)(i + shiftVelocity, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + shiftVelocity, jj) += MIn(ii, jj);
                (*M)(ii, j + shiftVelocity) += MOut(ii, jj);
            }
        }
    }

    // Computes the integral {u dot n}{v dot n}

    scalarVector_Type otherFractureEtaNormal ( fracture->getMeshFEMPressure().nb_dof(), 0. );

    for ( size_type i = 0; i < fracture->getMeshFEMPressure().nb_dof(); ++i )
    {
        base_node nodo = fracture->getMeshFEMPressure().point_of_basic_dof( i );
        otherFractureEtaNormal [ i ] = otherFracture->getData().getEtaNormal () *
                    otherFracture->getData().etaNormalDistribution ( nodo ) / fracture->getData().getThickness();
    }

    assemGam.set("w=data(#2);"
        "a=comp(vBase(#1).NonLin(#3).vBase(#1).NonLin(#3).Base(#2));"
        "M(#1,#1)+=a(:,j,j,:,i,i,k).w(k);");

    level_set_unit_normal nterm ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem(),
                                  fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->values() );

    getfem::mesh_im_level_set meshImLevel ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                            getfem::mesh_im_level_set::INTEGRATE_BOUNDARY );

    meshImLevel.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevel.set_simplex_im ( intTypeIM );

    // Assign the mesh integration method
    assemGam.push_mi ( meshImLevel );

    // Assign the mesh finite element space
    assemGam.push_mf ( fracture->getMeshFEMVelocity() );
    assemGam.push_mf ( fracture->getMeshFEMPressure() );

    // Assign the non linear term
    assemGam.push_nonlinear_term ( &nterm );

    // Assign the mesh finite element space for the coefficients
    assemGam.push_mf ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem() );

    // Assign the coefficients
    assemGam.push_data ( otherFractureEtaNormal );

    // Set the matrices to save the evaluations
    assemGam.push_mat ( Gamma  );

    // Computes the matrices
    assemGam.assembly ( cutRegionFlag );

    // Add the extended degrees of freedom
    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
/*        const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
        base_node pointMapped(0,0);
        pointMapped[0] = pointFlat[0];
        pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
        const scalar_type levelSetValue1 = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );
*/
        for ( size_type j = 0; j < extendedNumVelocity; ++j )
        {
            const size_type jj = extendedVelocity [ j ];
/*            const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(jj);
            base_node pointMapped(0,0);
            pointMapped[0] = pointFlat[0];
            pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
            const scalar_type levelSetValue2 = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );
*/
            (*M)(ii, jj) += 0.25 * Gamma(ii, jj);
            (*M)(i + shiftVelocity, j + shiftVelocity) += 0.25 * Gamma(ii, jj);
            (*M)(i + shiftVelocity, jj) += 0.25 * Gamma(ii, jj);
            (*M)(ii, j + shiftVelocity) += 0.25 * Gamma(ii, jj);
        }
    }

} // darcy_A11F

void divHdiv ( sparseMatrixPtr_Type& M,
               const MeshHandlerPtr_Type& mesh,
               const size_type& uncutRegionFlag )
{

    const size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();
    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize(M_, shiftVelocity, shiftPressure);
    gmm::clear(M_);

    // Volume integration
    getfem::generic_assembly assem;

    assem.set("M$1(#1,#2)+=-comp(vGrad(#1).Base(#2))"
        "(:,i,i, :);");

    // Assign the M_mediumMesh integration method
    assem.push_mi(mesh->getIntegrationMethodVector());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(mesh->getMeshFEMVector());
    assem.push_mf(mesh->getMeshFEMScalar());

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly(mesh->getRegion(uncutRegionFlag));

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftPressure; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    }

    cout << "DARCY :: operator b(volumic)     [OK]" << endl;

} // divHdiv

void divHdiv ( sparseMatrixPtr_Type& M,
               const MeshHandlerPtr_Type& mesh,
               const FractureHandlerPtr_Type& fracture,
               const size_type& cutRegionFlag )
{
    const scalar_type fractureID = fracture->getId();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const sizeVector_Type& extendedDOFScalar = mesh->getExtendedDOFScalar(
            fractureID);
    const size_type shiftExtendedScalar = extendedDOFScalar.size();
    const sizeVector_Type& extendedDOFVector = mesh->getExtendedDOFVector(
            fractureID);
    const size_type shiftExtendedVector = extendedDOFVector.size();

    size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();
    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();

    sparseMatrix_Type MIn, MOut;
    gmm::resize(MIn, shiftVelocity, shiftPressure);
    gmm::clear(MIn);
    gmm::resize(MOut, shiftVelocity, shiftPressure);
    gmm::clear(MOut);

    shiftPressure += mesh->getCountExtendedDOFScalar(fractureID - 1);
    shiftVelocity += mesh->getCountExtendedDOFVector(fractureID - 1);

    // Volume integration
    getfem::generic_assembly assemIn, assemOut;

    assemIn.set("M$1(#1,#2)+=-comp(vGrad(#1).Base(#2))"
        "(:,i,i, :);");

    assemOut.set("M$1(#1,#2)+=-comp(vGrad(#1).Base(#2))"
        "(:,i,i, :);");

    // Assign the mesh integration method
    assemIn.push_mi(levelSet->getIntegrationMethodInside());
    assemOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the mesh finite element space
    assemIn.push_mf(mesh->getMeshFEMVector());
    assemOut.push_mf(mesh->getMeshFEMVector());
    assemIn.push_mf(mesh->getMeshFEMScalar());
    assemOut.push_mf(mesh->getMeshFEMScalar());

    // Set the matrices to save the evaluations
    assemIn.push_mat(MIn);
    assemOut.push_mat(MOut);

    // Computes the matrices
    assemIn.assembly(mesh->getRegion(cutRegionFlag));
    assemOut.assembly(mesh->getRegion(cutRegionFlag));

    // Update the extended degrees of freedom
    for ( size_type i = 0; i < shiftExtendedVector; ++i )
    {
        size_type ii = extendedDOFVector [ i ];

        for ( size_type j = 0; j < shiftExtendedScalar; ++j )
        {
            size_type jj = extendedDOFScalar [ j ];

            if ( (levelSet->getDOFValue(ii) < 0)
                    && (levelSet->getBaricenterValue(jj) < 0) )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn(ii, jj);
                (*M)(i + shiftVelocity, j + shiftPressure) += MOut(ii, jj);
            }
            else if ( (levelSet->getDOFValue(ii) >= 0)
                    && (levelSet->getBaricenterValue(jj) >= 0) )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + shiftVelocity, j + shiftPressure) += MIn(ii, jj);
            }
            else if ( (levelSet->getDOFValue(ii) < 0)
                    && (levelSet->getBaricenterValue(jj) >= 0) )
            {
                // i is In, j is Out
                (*M)(ii, j + shiftPressure) += MIn(ii, jj);
                (*M)(i + shiftVelocity, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + shiftVelocity, jj) += MIn(ii, jj);
                (*M)(ii, j + shiftPressure) += MOut(ii, jj);
            }
        }

    }

    cout << "DARCY :: operator b(volumic)     [OK]" << endl;

} // divHdiv

void divHdiv ( sparseMatrixPtr_Type& M, const IntersectData_Type& intersect, const MeshHandlerPtr_Type& mesh )
{
    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();
    sparseMatrix_Type Matrix;
    gmm::resize(Matrix, shiftVelocity, shiftPressure );

    // construct the integration method based on the meshLevelSet
    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_TRIANGLE(6)" );
    getfem::mesh_im_level_set meshImLevelSet ( mesh->getMeshLevelSet(), getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );
    meshImLevelSet.set_integration_method ( mesh->getMesh().convex_index(), intTypeIM );
    meshImLevelSet.set_simplex_im ( intTypeIM );

    getfem::mesh_region meshElement;
    meshElement.add ( intersect.getElementID() );
    sizeVector_Type levelSetID ( intersect.getNumFractures(), 0 );
    for ( size_type f = 0; f < levelSetID.size(); ++f )
    {
        levelSetID [f] = intersect.getFracture(f)->getId();
    }

    getfem::mesh_fem::ind_dof_ct indexDOFVel = mesh->getMeshFEMVector().ind_basic_dof_of_element ( intersect.getElementID() );
    size_type indexDOFScal = mesh->getMeshFEMScalar().ind_basic_dof_of_element ( intersect.getElementID() )[0];
    for ( size_type regionID = 0; regionID < intersect.getRegionActive().size(); ++regionID )
    {
        gmm::clear ( Matrix );
        generic_assembly assembly;
        assembly.set("M$1(#1,#2)+=-comp(vGrad(#1).Base(#2))" "(:,i,i, :);");

        // Set the operation bewteen the level sets
        const std::string activeRegion = intersect.getRegionActive()[regionID];

        const std::string operation = getOperation ( activeRegion, levelSetID );

        meshImLevelSet.set_level_set_boolean_operations ( operation );
        meshImLevelSet.adapt();

        assembly.push_mi ( meshImLevelSet );
        assembly.push_mf ( mesh->getMeshFEMVector() );
        assembly.push_mf ( mesh->getMeshFEMScalar() );
        assembly.push_mat ( Matrix );
        assembly.assembly ( meshElement );

        // save the solution in M
        const size_type shiftJ = intersect.getDOFPressure ( regionID );

        // loop on rows
        for ( size_type i = 0; i < indexDOFVel.size(); ++i )
        {
            const size_type shiftI = intersect.getDOFVelocity ( i, regionID );

            // loop on column
            (*M)( shiftI, shiftJ ) += Matrix ( indexDOFVel[i], indexDOFScal );
        }
    }

} // divHdiv

//lo stesso per la frattura
void darcy_A12F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const size_type& uncutRegionFlag )
{
    const size_type velocityShift = fracture->getMeshFEMVelocity().nb_dof();
    const size_type pressureShift = fracture->getMeshFEMPressure().nb_dof();

    sparseMatrix_Type M_;

    gmm::resize(M_, velocityShift, pressureShift);
    gmm::clear(M_);

    // Volume integration
    getfem::generic_assembly assem;

    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
            "(:, i,i,:);");
    }
    else
    {
        assem.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
            "(:,i,i, :);");
    }

    // Assign the M_mediumMesh integration method
    assem.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(fracture->getMeshFEMVelocity());
    assem.push_mf(fracture->getMeshFEMPressure());

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly ( uncutRegionFlag );

    for ( size_type i = 0; i < velocityShift; ++i )
    {
        for ( size_type j = 0; j < pressureShift; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operator b(volumic)     [OK]" << endl;

} // darcy_A12F

//lo stesso per la frattura
void darcy_A12F ( sparseMatrixPtr_Type& M,
                  const FractureHandlerPtr_Type& fracture,
                  const FractureHandlerPtr_Type& otherFracture,
                  const size_type& cutRegionFlag )
{
    const size_type velocityShift = fracture->getMeshFEMVelocity().nb_dof();
    const size_type pressureShift = fracture->getMeshFEMPressure().nb_dof();

    sparseMatrix_Type MIn, MOut;

    gmm::resize(MIn, velocityShift, pressureShift);
    gmm::clear(MIn);
    gmm::resize(MOut, velocityShift, pressureShift);
    gmm::clear(MOut);

    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    // Volume integration
    getfem::generic_assembly assemIn, assemOut;

    assemIn.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
                "(:, i,i,:);");
    assemOut.set("M(#1,#2)+=-comp(vGrad(#1).Base(#2))"
                 "(:, i,i,:);");

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevelSetOut ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                  getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );

    getfem::mesh_im_level_set meshImLevelSetIn ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                 getfem::mesh_im_level_set::INTEGRATE_INSIDE );

    meshImLevelSetOut.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );
    meshImLevelSetIn.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevelSetOut.set_simplex_im ( intTypeIM );
    meshImLevelSetIn.set_simplex_im ( intTypeIM );

    // Assign the mesh integration method
    assemIn.push_mi ( meshImLevelSetIn );
    assemOut.push_mi ( meshImLevelSetOut );

    // Assign the M_mediumMesh finite element space
    assemIn.push_mf ( fracture->getMeshFEMVelocity() );
    assemIn.push_mf ( fracture->getMeshFEMPressure() );

    assemOut.push_mf ( fracture->getMeshFEMVelocity() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Set the matrices to save the evaluations
    assemIn.push_mat ( MIn );
    assemOut.push_mat ( MOut );

    // Computes the matrices
    assemIn.assembly ( cutRegionFlag );
    assemOut.assembly ( cutRegionFlag );

    // Update the extended degrees of freedom
    const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();
    const size_type extendedNumVelocity = fracture->getNumExtendedVelocity();
    const sizeVector_Type& extendedPressure = fracture->getExtendedPressure();
    const size_type extendedNumPressure = fracture->getNumExtendedPressure();

    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
        const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
        base_node pointMapped(0,0);
        pointMapped[0] = pointFlat[0];
        pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
        const scalar_type levelSetValue1 = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );

        for ( size_type j = 0; j < extendedNumPressure; ++j )
        {
            const size_type jj = extendedPressure [ j ];
            const base_node pointFlat = fracture->getMeshFEMPressure().point_of_basic_dof(jj);
            base_node pointMapped(0,0);
            pointMapped[0] = pointFlat[0];
            pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
            const scalar_type levelSetValue2 = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );

            if ( levelSetValue1 < 0 && levelSetValue2 < 0 )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn ( ii, jj );
                (*M)(i + velocityShift, j + pressureShift) += MOut ( ii, jj );
            }

            else if ( levelSetValue1 >= 0 && levelSetValue2 >= 0 )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + velocityShift, j + pressureShift) += MIn(ii, jj);
            }
            else if ( levelSetValue1 < 0 && levelSetValue2 >= 0 )
            {
                // i is In, j is Out
                (*M)(ii, j + pressureShift) += MIn(ii, jj);
                (*M)(i + velocityShift, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + velocityShift, jj) += MIn(ii, jj);
                (*M)(ii, j + pressureShift) += MOut(ii, jj);
            }
        }
    }

} // darcy_A12F

//termine di trasporto
void advectHdiv ( sparseMatrixPtr_Type& M,
                  const MeshHandlerPtr_Type& mesh,
                  const scalarVectorContainer_Type& mediumMu,
                  const scalarVectorPtr_Type& velocity,
                  const size_type& uncutRegionFlag )
{
    const size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();
    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();

    sparseMatrix_Type M_;
    gmm::resize(M_, shiftVelocity, shiftPressure);
    gmm::clear(M_);

    // Volume integration
    getfem::generic_assembly assem;

    /// ERRORE SE USO K TENSORIALE
    assem.set("mediumMu=data$1(#1);" "velocity=data$2(#2);"
        "M(#2,#1)+=-comp(Base(#1).vBase(#2).vBase(#2).Base(#1))"
        "(j,l,i,:,i,:).mediumMu(j).velocity(l);");

    // Assign the M_mediumMesh integration method
    assem.push_mi(mesh->getIntegrationMethodVector());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(mesh->getMeshFEMScalar());
    assem.push_mf(mesh->getMeshFEMVector());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem.push_mf(mesh->getMeshFEMScalar());
    assem.push_mf(mesh->getMeshFEMVector());

    // Assign the coefficients
    assem.push_data(mediumMu[0]);
    assem.push_data(*velocity);

    // Set the matrices to save the evaluations
    assem.push_mat(M_);

    // Computes the matrices
    assem.assembly(mesh->getRegion(uncutRegionFlag));

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        for ( size_type j = 0; j < shiftPressure; ++j )
        {
            (*M)(i, j) = M_(i, j);
        }
    };

    cout << "DARCY :: operatore di trasporto [OK]" << endl;

} // darcy_D

//termine di trasporto
void advectHdiv ( sparseMatrixPtr_Type& M,
                  const MeshHandlerPtr_Type& mesh,
                  const FractureHandlerPtr_Type& fracture,
                  const scalarVectorContainer_Type& mediumMu,
                  const scalarVectorPtr_Type& velocityIn,
                  const scalarVectorPtr_Type& velocityOut,
                  const size_type& cutRegionFlag )
{
    const scalar_type fractureID = fracture->getId();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();
    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    const sizeVector_Type& extendedDOFScalar = mesh->getExtendedDOFScalar(
            fractureID);
    const size_type shiftExtendedScalar = extendedDOFScalar.size();

    const sizeVector_Type& extendedDOFVector = mesh->getExtendedDOFVector(
            fractureID);
    const size_type shiftExtendedVector = extendedDOFVector.size();

    sparseMatrix_Type MIn, MOut;
    gmm::resize(MIn, shiftVelocity, shiftPressure);
    gmm::clear(MIn);
    gmm::resize(MOut, shiftVelocity, shiftPressure);
    gmm::clear(MOut);

    shiftPressure += mesh->getCountExtendedDOFScalar(fractureID - 1);
    shiftVelocity += mesh->getCountExtendedDOFVector(fractureID - 1);

    // Volume integration
    getfem::generic_assembly assemIn, assemOut;


    // ERRORE SE USO K TENSORE
    assemIn.set("mediumMu=data$1(#1);" "velocity=data$2(#2);"
        "M(#2,#1)+=-comp(Base(#1).vBase(#2).vBase(#2).Base(#1))"
        "(j,l,i,:,i,:).mediumMu(j).velocity(l);");

    assemOut.set("mediumMu=data$1(#1);" "velocity=data$2(#2);"
        "M(#2,#1)+=-comp(Base(#1).vBase(#2).vBase(#2).Base(#1))"
        "(j,l,i,:,i,:).mediumMu(j).velocity(l);");

    // Assign the M_mediumMesh integration method
    assemIn.push_mi(levelSet->getIntegrationMethodInside());
    assemOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the M_mediumMesh finite element space
    assemIn.push_mf(mesh->getMeshFEMScalar());
    assemOut.push_mf(mesh->getMeshFEMScalar());
    assemIn.push_mf(mesh->getMeshFEMVector());
    assemOut.push_mf(mesh->getMeshFEMVector());

    // Assign the M_mediumMesh finite element space for the coefficients
    assemIn.push_mf(mesh->getMeshFEMScalar());
    assemOut.push_mf(mesh->getMeshFEMScalar());
    assemIn.push_mf(mesh->getMeshFEMVector());
    assemOut.push_mf(mesh->getMeshFEMVector());

    // Assign the coefficients
    assemIn.push_data(mediumMu[0]);
    assemOut.push_data(mediumMu[0]);
    assemIn.push_data(*velocityIn);
    assemOut.push_data(*velocityOut);

    // Set the matrices to save the evaluations
    assemIn.push_mat(MIn);
    assemOut.push_mat(MOut);

    // Computes the matrices
    assemIn.assembly(mesh->getRegion(cutRegionFlag));
    assemOut.assembly(mesh->getRegion(cutRegionFlag));

    // Update the extended degrees of freedom
    for ( size_type i = 0; i < shiftExtendedVector; ++i )
    {
        size_type ii = extendedDOFVector [ i ];

        for ( size_type j = 0; j < shiftExtendedScalar; ++j )
        {
            size_type jj = extendedDOFScalar [ j ];
            if ( (levelSet->getDOFValue(ii) < 0)
                    && (levelSet->getBaricenterValue(jj) < 0) )
            {
                // i and j are both In
                (*M)(ii, jj) += MIn(ii, jj);
                (*M)(i + shiftVelocity, j + shiftPressure) += MOut(ii, jj);
            }
            else if ( (levelSet->getDOFValue(ii) >= 0)
                    && (levelSet->getBaricenterValue(jj) >= 0) )
            {
                // i and j are both Out
                (*M)(ii, jj) += MOut(ii, jj);
                (*M)(i + shiftVelocity, j + shiftPressure) += MIn(ii, jj);
            }
            else if ( (levelSet->getDOFValue(ii) < 0)
                    && (levelSet->getBaricenterValue(jj) >= 0) )
            {
                // i is In, j is Out
                (*M)(ii, j + shiftPressure) += MIn(ii, jj);
                (*M)(i + shiftVelocity, jj) += MOut(ii, jj);
            }
            else
            {
                // i is Out, j is In
                (*M)(i + shiftVelocity, jj) += MIn(ii, jj);
                (*M)(ii, j + shiftPressure) += MOut(ii, jj);
            }
        }

    }

    cout << "DARCY :: operatore di trasporto [OK]" << endl;

} // darcy_D

// termine di trasporto nella frattura
void darcy_DF ( sparseMatrixPtr_Type& M,
                const getfem::mesh_im& fractureIntegrationMethodVelocity,
                const getfem::mesh_fem& fractureMeshFEMVelocity,
                const getfem::mesh_fem& fractureMeshFEMPressure,
                const scalarVector_Type& fractureMuTangential,
                const scalarVector_Type& fractureMagnificationMapFactor1,
                const scalarVectorPtr_Type& velocity,
                const mesh_region& uncutRegion )
{
    scalarVector_Type invF(fractureMagnificationMapFactor1);

    // Volume integration
    getfem::generic_assembly assem1;
    assem1.set(
            "invF=data$1(#1);" "fractureMuTangential=data$2(#1);" "velocity=data$3(#2);"
                "M(#2,#1)+=-comp(Base(#1).Base(#1).vBase(#2).vBase(#2).Base(#1))"
                "(k,j,l,i,:,i,:).invF(k).fractureMuTangential(j).velocity(l);");

    // Assign the M_mediumMesh integration method
    assem1.push_mi(fractureIntegrationMethodVelocity);

    // Assign the M_mediumMesh finite element space
    assem1.push_mf(fractureMeshFEMPressure);
    assem1.push_mf(fractureMeshFEMVelocity);

    // Assign the M_mediumMesh finite element space for the coefficients
    assem1.push_mf(fractureMeshFEMPressure);
    assem1.push_mf(fractureMeshFEMPressure);
    assem1.push_mf(fractureMeshFEMVelocity);

    // Assign the coefficients
    assem1.push_data(invF);
    assem1.push_data(fractureMuTangential);
    assem1.push_data(*velocity);

    // Set the matrices to save the evaluations
    assem1.push_mat(*M);

    // Computes the matrices
    assem1.assembly(uncutRegion);

    cout << "DARCY :: operatore di trasporto [OK]" << endl;

} // darcy_DF

//questa fa la parte di termine noto che contiene le condizioni al contorno
//il termine di pressione entra come condizione naturale, quello di velocità con una penalizzazione di Nitzche
void stressRHS ( scalarVectorPtr_Type &Bstress,
                 const MeshHandlerPtr_Type& mesh,
                 const BCHandlerPtr_Type& bcHandler,
                 const scalarVectorPtr_Type& Pneumann )
{

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();

    // ----------------- Ext Stress - dato di pressione ---------------

    scalarVector_Type Bs(shiftVelocity, 0.);

    //assemblo pd tau dot n
    getfem::generic_assembly assemb;

    assemb.set("p=data(#2);"
        "V(#1)+=-comp(vBase(#1).Normal().Base(#2))"
        "(:,k, k, h).p(h);");

    // Assign the M_mediumMesh integration method
    assemb.push_mi(mesh->getIntegrationMethodVector());

    // Assign the M_mediumMesh finite element space
    assemb.push_mf(mesh->getMeshFEMVector());

    // Assign the M_mediumMesh finite element space for the coefficients
    assemb.push_mf(bcHandler->getMediumBC()->getMeshFEM());

    // Assign the coefficients
    assemb.push_data(*Pneumann);

    // Set the vectors to save the evaluations
    assemb.push_vec(Bs);

    // Assemble in each sub region
    const size_type shiftNeumann = bcHandler->getNeumannUncut().size();
    for ( size_type bndID = 0; bndID < shiftNeumann; bndID++ )
    {
        const size_type val = bcHandler->getNeumannUncut(bndID);
        assemb.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(val));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bstress) [ i ] = Bs [ i ];
    }

} // darcy_data

//questa fa la parte di termine noto che contiene le condizioni al contorno
//il termine di pressione entra come condizione naturale, quello di velocità con una penalizzazione di Nitzche
void nitscheRHS ( scalarVectorPtr_Type& Bvel,
                  const MeshHandlerPtr_Type& mesh,
                  const BCHandlerPtr_Type& bcHandler,
                  const scalar_type& eta,
                  const scalar_type& gammaU,
                  const scalarVectorPtr_Type& v_diri )
{

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    const size_type shiftCoefficients = mesh->getMeshFEMCoefficients().nb_dof();

    // ----------------- Penalty    ---------------

    scalarVector_Type etaGammaUinvh(shiftCoefficients);
    sizeVector_Type Bvel_tot(shiftVelocity, 0.);
    sizeVector_Type BvelIn(shiftVelocity, 0.);
    sizeVector_Type BvelOut(shiftVelocity, 0.);

    for ( size_type i = 0; i < shiftCoefficients; ++i )
    {
        etaGammaUinvh [ i ] = eta * gammaU * mesh->getInverseMeshSizeDOF(i);
    }

    getfem::generic_assembly assem;

    //questo è gamma v dot n tau dot n.
    //come dato di dirichlet nota bene conosco direttamente - dal main - la v NORMALE

    assem.set("gamma=data$1(#2);" "vel=data$2(#2);"
        "t=comp(Base(#2).vBase(#1).Normal().Base(#2));"
        "V(#1)+=(t(m, :,j, j, k).gamma(k).vel(m));");

    // Assign the M_mediumMesh integration method
    assem.push_mi(mesh->getIntegrationMethodVector());

    // Assign the M_mediumMesh finite element space
    assem.push_mf(mesh->getMeshFEMVector());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem.push_mf(mesh->getMeshFEMCoefficients());
    assem.push_mf(mesh->getMeshFEMCoefficients());

    scalarVector_Type eta_provv_Un(etaGammaUinvh);

    // Assign the coefficients
    assem.push_data(eta_provv_Un);
    assem.push_data(*v_diri);

    // Set the vector to save the evaluations
    assem.push_vec(Bvel_tot);

    // Assemble in each sub region
    const size_type shiftDirichlet = bcHandler->getDirichletUncut().size();
    for ( size_type bndID = 0; bndID < shiftDirichlet; bndID++ )
    {
        const size_type val = bcHandler->getDirichletUncut(bndID);
        assem.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(val));

    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bvel) [ i ] += Bvel_tot [ i ];
    }

    cout << "DARCY :: DATA (penal. bound.)    [OK]" << endl;

} // darcy_data

//questa fa la parte di termine noto che contiene le condizioni al contorno
//il termine di pressione entra come condizione naturale, quello di velocità con una penalizzazione di Nitzche
void stressRHS ( scalarVectorPtr_Type &Bstress,
                 const MeshHandlerPtr_Type& mesh,
                 const BCHandlerPtr_Type& bcHandler,
                 const FractureHandlerPtr_Type& fracture,
                 const size_type& fractureShift,
                 const scalarVectorPtr_Type& PneumannIn,
                 const scalarVectorPtr_Type& PneumannOut )
{

    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const scalar_type fractureID = fracture->getId();
    const sizeVector_Type& extendedDOF = mesh->getExtendedDOFVector(fractureID);
    const size_type shiftExtended = extendedDOF.size();
    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();

    // ----------------- Ext Stress - dato di pressione ---------------

    scalarVector_Type BsIn(shiftVelocity, 0.), BsOut(shiftVelocity, 0.);

    shiftVelocity += fractureShift;

    //assemblo pd tau dot n
    getfem::generic_assembly assembIn, assembOut;

    assembIn.set("p=data(#2);"
        "V(#1)+=-comp(vBase(#1).Normal().Base(#2))"
        "(:,k, k, h).p(h);");

    assembOut.set("p=data(#2);"
        "V(#1)+=-comp(vBase(#1).Normal().Base(#2))"
        "(:,k, k, h).p(h);");

    // Assign the mesh integration method
    assembIn.push_mi(levelSet->getIntegrationMethodInside());
    assembOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the mesh finite element space
    assembIn.push_mf(mesh->getMeshFEMVector());
    assembOut.push_mf(mesh->getMeshFEMVector());

    // Assign the mesh finite element space for the coefficients
    assembIn.push_mf(bcHandler->getMediumBC()->getMeshFEM());
    assembOut.push_mf(bcHandler->getMediumBC()->getMeshFEM());

    // Assign the coefficients
    assembIn.push_data(*PneumannIn);
    assembOut.push_data(*PneumannOut);

    // Set the vectors to save the evaluations
    assembIn.push_vec(BsIn);
    assembOut.push_vec(BsOut);

    // Assemble in each sub region
    const size_type shiftNeumann = bcHandler->getNeumannCut(fractureID).size();
    for ( size_type bndID = 0; bndID < shiftNeumann; bndID++ )
    {
        const size_type val = bcHandler->getNeumannCut(fractureID, bndID);
        assembIn.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(val));
        assembOut.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(val));
    }

    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = extendedDOF [ i ];

        if ( levelSet->getDOFValue(ii) < 0 )
        {
            (*Bstress) [ ii ] += BsIn [ ii ];
            (*Bstress) [ i + shiftVelocity ] += BsOut [ ii ];
        }
        else
        {
            (*Bstress) [ ii ] += BsOut [ ii ];
            (*Bstress) [ i + shiftVelocity ] += BsIn [ ii ];
        }
    }

} // darcy_data

//questa fa la parte di termine noto che contiene le condizioni al contorno
//il termine di pressione entra come condizione naturale, quello di velocità con una penalizzazione di Nitzche
void nitscheRHS ( scalarVectorPtr_Type &Bvel,
                  const MeshHandlerPtr_Type& mesh,
                  const BCHandlerPtr_Type& bcHandler,
                  const FractureHandlerPtr_Type& fracture,
                  const size_type& fractureShift,
                  const scalar_type& invK,
                  const scalar_type& gammaU,
                  const scalarVectorPtr_Type& v_diri )
{
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const scalar_type fractureID = fracture->getId();
    const sizeVector_Type& extendedDOF = mesh->getExtendedDOFVector(fractureID);
    const size_type shiftExtended = extendedDOF.size();

    size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    const size_type shiftCoefficients = mesh->getMeshFEMCoefficients().nb_dof();

    // ----------------- Penalty    ---------------

    scalarVector_Type etaGammaUinvh(shiftCoefficients);
    sizeVector_Type BvelIn(shiftVelocity, 0.), BvelOut(shiftVelocity, 0.);

    shiftVelocity += fractureShift;

    for ( size_type i = 0; i < shiftCoefficients; ++i )
    {
        etaGammaUinvh [ i ] = invK * gammaU * mesh->getInverseMeshSizeDOF(i);
    }

    getfem::generic_assembly assemIn, assemOut;

    //questo è gamma v dot n tau dot n.
    //come dato di dirichlet nota bene conosco direttamente - dal main - la v NORMALE
    assemIn.set("gamma=data$1(#2);" "vel=data$2(#2);"
        "t=comp(Base(#2).vBase(#1).Normal().Base(#2));"
        "V(#1)+=(t(m, :,j, j, k).gamma(k).vel(m));");

    assemOut.set("gamma=data$1(#2);" "vel=data$2(#2);"
        "t=comp(Base(#2).vBase(#1).Normal().Base(#2));"
        "V(#1)+=(t(m, :,j, j, k).gamma(k).vel(m));");

    // Assign the mesh integration method
    assemIn.push_mi(levelSet->getIntegrationMethodInside());
    assemOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the mesh finite element space
    assemIn.push_mf(mesh->getMeshFEMVector());
    assemOut.push_mf(mesh->getMeshFEMVector());

    // Assign the mesh finite element space for the coefficients
    assemIn.push_mf(mesh->getMeshFEMCoefficients());
    assemOut.push_mf(mesh->getMeshFEMCoefficients());
    assemIn.push_mf(mesh->getMeshFEMCoefficients());
    assemOut.push_mf(mesh->getMeshFEMCoefficients());

    scalarVector_Type eta_provv_Cut(etaGammaUinvh);

    // Assign the coefficients
    assemIn.push_data(eta_provv_Cut);
    assemOut.push_data(eta_provv_Cut);
    assemIn.push_data(*v_diri);
    assemOut.push_data(*v_diri);

    // Set the vector to save the evaluations
    assemIn.push_vec(BvelIn);
    assemOut.push_vec(BvelOut);

    // Assemble in each sub region
    const size_type shiftDirichlet =
            bcHandler->getDirichletCut(fractureID).size();
    for ( size_type bndID = 0; bndID < shiftDirichlet; bndID++ )
    {
        const size_type val = bcHandler->getDirichletCut(fractureID, bndID);
        assemIn.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(val));
        assemOut.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(val));
    }

    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = extendedDOF [ i ];

        if ( levelSet->getDOFValue(ii) < 0 )
        {
            (*Bvel) [ ii ] += BvelIn [ ii ];
            (*Bvel) [ i + shiftVelocity ] += BvelOut [ ii ];
        }
        else
        {
            (*Bvel) [ ii ] += BvelOut [ ii ];
            (*Bvel) [ i + shiftVelocity ] += BvelIn [ ii ];
        }
    }
    cout << "DARCY :: DATA (penal. bound.)    [OK]" << endl;

} // darcy_data

//stesso lavoro con la frattura - più semplice
void darcy_dataF ( scalarVectorPtr_Type &Bstress,
                   scalarVectorPtr_Type &Bvel,
                   const BCHandlerPtr_Type& bcHandler,
                   const FractureHandlerPtr_Type& fracture,
                   const scalar_type& gammaU,
                   const scalar_type& invK,
                   const scalarVectorPtr_Type& Pneumann,
                   const scalarVectorPtr_Type& v_diri )
{

    // ----------------- Penalty    ---------------

    const scalar_type fractureID = fracture->getId();
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    const size_type shiftCoefficinents =
            fracture->getMeshFEMPressure().nb_dof();

    scalarVector_Type etaGammaUinvh(shiftCoefficinents, 0.), Bvel_tot(
            shiftVelocity, 0.);

    for ( size_type i = 0; i < shiftCoefficinents; ++i )
    {
        etaGammaUinvh [ i ] = invK * gammaU * fracture->getInverseMeshSize(i);
    }

    getfem::generic_assembly assem2;

    assem2.set("gamma=data$1(#2);" "vel=data$2(#2);"
        "t=comp(Base(#2).vBase(#1).Normal().Base(#2));"
        "V(#1)+=(t(m, :,j, j, k).gamma(k).vel(m));");

    // Assign the M_mediumMesh integration method
    assem2.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assem2.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem2.push_mf(fracture->getMeshFEMPressure());
    assem2.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem2.push_data(etaGammaUinvh);
    assem2.push_data(*v_diri);

    // Set the matrices to save the evaluations
    assem2.push_vec(Bvel_tot);

    // Assemble in each sub region
    const size_type shiftDirichlet =
            bcHandler->getFractureBC(fractureID)->getDirichlet().size();
    for ( size_type bndID = 0; bndID < shiftDirichlet; bndID++ )
    {
        const size_type val =
                bcHandler->getFractureBC(fractureID)->getDirichlet(bndID);
        assem2.assembly(
                fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region(
                        val));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bvel) [ i ] += Bvel_tot [ i ];
    }

    cout << "DARCY :: DATA (penal. bound.)    [OK]" << endl;

    // ----------------- Ext Stress ---------------

    const size_type shiftMapFactor1 =
            fracture->getMagnificationMapFactor1().size();
    const size_type shiftMapFactor2 =
            fracture->getMagnificationMapFactor2().size();

    scalarVector_Type Bs(shiftVelocity, 0.), coefx(shiftMapFactor1);
    scalarVector_Type coefy(shiftMapFactor2);

    getfem::generic_assembly assemb;

    if ( fracture->getMeshFEMVelocity().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor1; ++i )
        {
            coefx [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
        }
        assemb.set("p=data$1(#2);"
            "V(#1)+=-comp(vBase(#1).Base(#2))"
            "(:,1, h).p(h)");
    }
    else
    {

        //attenzione, quando integro sulla frattura devo tenere conto che quella vera non è flat quindi le lunghezze/aree devono essere convertite
        for ( size_type i = 0; i < shiftMapFactor1; ++i )
        {
            coefx [ i ] = 1 / fracture->getMagnificationMapFactor1(i);
            coefy [ i ] = 1 / fracture->getMagnificationMapFactor2(i);
        }
        assemb.set("p=data$1(#2);"
            "V(#1)+=-comp(vBase(#1).Normal().Base(#2))"
            "(:,k, k, h).p(h)");
    }

    // Assign the M_mediumMesh integration method
    assemb.push_mi(fracture->getIntegrationMethodVelocity());

    // Assign the M_mediumMesh finite element space
    assemb.push_mf(fracture->getMeshFEMVelocity());

    // Assign the M_mediumMesh finite element space for the coefficients
    assemb.push_mf(bcHandler->getFractureBC(fracture->getId())->getMeshFEM());

    // Assign the coefficients
    assemb.push_data(*Pneumann);

    // Set the vector to save the evaluations
    assemb.push_vec(Bs);

    // Assemble in each sub region
    const size_type shiftNeumann =
            bcHandler->getFractureBC(fractureID)->getNeumann().size();
    for ( size_type bndID = 0; bndID < shiftNeumann; bndID++ )
    {
        const size_type val = bcHandler->getFractureBC(fractureID)->getNeumann(
                bndID);
        assemb.assembly(
                fracture->getMeshFEMVelocity().linked_mesh().get_mpi_sub_region(
                        val));
    }

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        (*Bstress) [ i ] = Bs [ i ];
    }

} // darcy_dataF

//questo tiene conto del termine sorgenre
void sourceL2 ( scalarVectorPtr_Type& D,
                const scalarVectorPtr_Type& source,
                const MeshHandlerPtr_Type& mesh,
                const size_type& uncutRegionFlag )
{

    const size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();

    scalarVector_Type D_(shiftPressure, 0.);

    generic_assembly assem_Source;

    assem_Source.set("w=data(#2);"
        "a=comp(Base(#1).Base(#2));"
        "V(#1)+=a(:, k).w(k)");

    // Assign the M_mediumMesh integration method
    assem_Source.push_mi(mesh->getIntegrationMethodScalar());

    // Assign the M_mediumMesh finite element space
    assem_Source.push_mf(mesh->getMeshFEMScalar());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_Source.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the coefficients
    assem_Source.push_data(*source);

    // Set the vectors to save the evaluations
    assem_Source.push_vec(D_);

    // Computes the matrices
    assem_Source.assembly(mesh->getRegion(uncutRegionFlag));

    gmm::add(D_, gmm::sub_vector(*D, gmm::sub_interval(0, shiftPressure)));

} // assembling_Source_Boundary

//questo tiene conto del termine sorgenre
void sourceL2 ( scalarVectorPtr_Type& D,
                const scalarVectorPtr_Type& source,
                const MeshHandlerPtr_Type& mesh,
                const FractureHandlerPtr_Type& fracture,
                const size_type& cutRegionFlag )
{
    const LevelSetHandlerPtr_Type& levelSet = fracture->getLevelSet();
    const scalar_type fractureID = fracture->getId();
    const sizeVector_Type& extendedDOF = mesh->getExtendedDOFScalar(fractureID);
    const size_type shiftExtended = extendedDOF.size();
    size_type shiftPressure = mesh->getMeshFEMScalar().nb_dof();

    scalarVector_Type DOut(shiftPressure, 0.), DIn(shiftPressure, 0.);
    shiftPressure += mesh->getCountExtendedDOFScalar(fractureID - 1);

    generic_assembly assem_SourceIn, assem_SourceOut;

    assem_SourceIn.set("w=data(#2);"
        "a=comp(Base(#1).Base(#2));"
        "V(#1)+=a(:, k).w(k)");

    assem_SourceOut.set("w=data(#2);"
        "a=comp(Base(#1).Base(#2));"
        "V(#1)+=a(:, k).w(k)");

    // Assign the M_mediumMesh integration method
    assem_SourceIn.push_mi(levelSet->getIntegrationMethodInside());
    assem_SourceOut.push_mi(levelSet->getIntegrationMethodOutside());

    // Assign the M_mediumMesh finite element space
    assem_SourceIn.push_mf(mesh->getMeshFEMScalar());
    assem_SourceOut.push_mf(mesh->getMeshFEMScalar());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_SourceIn.push_mf(mesh->getMeshFEMCoefficients());
    assem_SourceOut.push_mf(mesh->getMeshFEMCoefficients());

    // Assign the coefficients
    assem_SourceIn.push_data(*source);
    assem_SourceOut.push_data(*source);

    // Set the vectors to save the evaluations
    assem_SourceIn.push_vec(DIn);
    assem_SourceOut.push_vec(DOut);

    // Computes the matrices
    assem_SourceIn.assembly(mesh->getRegion(cutRegionFlag));
    assem_SourceOut.assembly(mesh->getRegion(cutRegionFlag));

    for ( size_type i = 0; i < shiftExtended; ++i )
    {
        size_type ii = extendedDOF [ i ];

        if ( levelSet->getDOFValue(ii) < 0 )
        {
            (*D) [ ii ] += DIn [ ii ];
            (*D) [ i + shiftPressure ] += DOut [ ii ];
        }
        else
        {
            (*D) [ ii ] += DOut [ ii ];
            (*D) [ i + shiftPressure ] += DIn [ ii ];
        }
    }

} // assembling_Source_Boundary

void sourceL2 ( scalarVectorPtr_Type& V, const scalarVectorPtr_Type& source,
                const IntersectData_Type& intersect, const MeshHandlerPtr_Type& mesh )
{
    size_type shiftScalar = mesh->getMeshFEMScalar().nb_dof();
    scalarVector_Type Vector( shiftScalar, 0 );

    getfem::mesh_fem::ind_dof_ct indexDOF = mesh->getMeshFEMScalar().ind_basic_dof_of_element ( intersect.getElementID() );

    for ( size_type regionID = 0; regionID < intersect.getRegionActive().size(); ++regionID )
    {
        gmm::clear ( Vector );
        generic_assembly assembly;
        assembly.set ( "w=data(#2);"
                       "a=comp(Base(#1).Base(#2));"
                       "V(#1)+=a(:, k).w(k)");

        // construct the integration method based on the meshLevelSet
        const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_TRIANGLE(6)" );
        getfem::mesh_im_level_set meshImLevelSet ( mesh->getMeshLevelSet(), getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );
        meshImLevelSet.set_integration_method ( mesh->getMesh().convex_index(), intTypeIM );
        meshImLevelSet.set_simplex_im ( intTypeIM );

        getfem::mesh_region meshElement;
        meshElement.add ( intersect.getElementID() );

        sizeVector_Type levelSetID ( intersect.getNumFractures(), 0 );
        for ( size_type f = 0; f < levelSetID.size(); ++f )
        {
            levelSetID [f] = intersect.getFracture(f)->getId();
        }

        assembly.push_mi ( meshImLevelSet );
        assembly.push_mf ( mesh->getMeshFEMScalar() );

        // Assign the mesh finite element space for the coefficients
        assembly.push_mf(mesh->getMeshFEMCoefficients());

        // Assign the coefficients
        assembly.push_data ( *source );

        assembly.push_vec ( Vector );

        // Set the operation bewteen the level sets
        const std::string activeRegion = intersect.getRegionActive()[regionID];

        const std::string operation = getOperation ( activeRegion, levelSetID );

        meshImLevelSet.set_level_set_boolean_operations ( operation );

        meshImLevelSet.adapt();

        assembly.assembly ( meshElement );

        const size_type shiftI = intersect.getDOFPressure ( regionID );

        (*V)[ shiftI ] += Vector [ indexDOF[0] ];

    }
}

//termine sorgente per la frattura
void assembling_Source_BoundaryF ( scalarVectorPtr_Type& D,
                                   const scalarVectorPtr_Type& source,
                                   const FractureHandlerPtr_Type& fracture,
                                   const size_type& uncutRegionFlag )
{

    const size_type shiftPressure = fracture->getMeshFEMPressure().nb_dof();
    const size_type shiftMapFactor =
            fracture->getMagnificationMapFactor1().size();

    scalarVector_Type D_(shiftPressure, 0.0), invF(shiftMapFactor, 0);

    generic_assembly assem_Source, assem_Vx, assem_Vy;

    if ( fracture->getMeshFEMPressure().get_qdim() == 1 )
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i));
        }

        assem_Source.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#2).Base(#2));"
            "V(#1)+=a(:, k,j).w(k).q(j)");

    }
    else
    {
        for ( size_type i = 0; i < shiftMapFactor; ++i )
        {
            invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i)
                    * fracture->getMagnificationMapFactor2(i));
        }

        assem_Source.set("w=data$1(#2);" "q=data$2(#2);"
            "a=comp(Base(#1).Base(#2).Base(#2));"
            "V(#1)+=a(:, k,j).w(k).q(j)");

    }

    // Assign the M_mediumMesh integration method
    assem_Source.push_mi(fracture->getIntegrationMethodPressure());

    // Assign the M_mediumMesh finite element space
    assem_Source.push_mf(fracture->getMeshFEMPressure());

    // Assign the M_mediumMesh finite element space for the coefficients
    assem_Source.push_mf(fracture->getMeshFEMPressure());
    assem_Source.push_mf(fracture->getMeshFEMPressure());

    // Assign the coefficients
    assem_Source.push_data(*source);
    assem_Source.push_data(invF);

    // Set the vector to save the evaluations
    assem_Source.push_vec(D_);

    // Computes the matrices
    assem_Source.assembly ( uncutRegionFlag );

    gmm::add(D_, gmm::sub_vector(*D, gmm::sub_interval(0, shiftPressure)));

} // assembling_Source_BoundaryF

//termine sorgente per la frattura
void assembling_SourceF ( scalarVectorPtr_Type& D,
                          const scalarVectorPtr_Type& source,
                          const FractureHandlerPtr_Type& fracture,
                          const FractureHandlerPtr_Type& otherFracture,
                          const size_type& cutRegionFlag )
{
    const size_type shiftPressure = fracture->getMeshFEMPressure().nb_dof();
    const size_type shiftMapFactor = fracture->getMagnificationMapFactor1().size();
    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    scalarVector_Type DIn ( shiftPressure, 0.0 ), DOut ( shiftPressure, 0.0 ),
                      invF ( shiftMapFactor, 0 );

    generic_assembly assemIn, assemOut;

    for ( size_type i = 0; i < shiftMapFactor; ++i )
    {
        invF [ i ] = 1.0 / (fracture->getMagnificationMapFactor1(i));
    }

    assemIn.set( "w=data$1(#2);" "q=data$2(#2);"
                 "a=comp(Base(#1).Base(#2).Base(#2));"
                 "V(#1)+=a(:, k,j).w(k).q(j)" );

    assemOut.set( "w=data$1(#2);" "q=data$2(#2);"
                  "a=comp(Base(#1).Base(#2).Base(#2));"
                  "V(#1)+=a(:, k,j).w(k).q(j)" );

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevelSetOut ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                  getfem::mesh_im_level_set::INTEGRATE_OUTSIDE );

    getfem::mesh_im_level_set meshImLevelSetIn ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                                 getfem::mesh_im_level_set::INTEGRATE_INSIDE );

    meshImLevelSetOut.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );
    meshImLevelSetIn.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevelSetOut.set_simplex_im ( intTypeIM );
    meshImLevelSetIn.set_simplex_im ( intTypeIM );

    // Assign the M_mediumMesh integration method
    assemIn.push_mi ( meshImLevelSetIn );
    assemOut.push_mi ( meshImLevelSetOut );

    // Assign the M_mediumMesh finite element space
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Assign the M_mediumMesh finite element space for the coefficients
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemIn.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );
    assemOut.push_mf ( fracture->getMeshFEMPressure() );

    // Assign the coefficients
    assemIn.push_data ( *source );
    assemIn.push_data ( invF );
    assemOut.push_data ( *source );
    assemOut.push_data ( invF );

    // Set the vector to save the evaluations
    assemIn.push_vec ( DIn );
    assemOut.push_vec ( DOut );

    // Computes the matrices
    assemIn.assembly ( cutRegionFlag );
    assemOut.assembly ( cutRegionFlag );

    const sizeVector_Type& extendedPressure = fracture->getExtendedPressure();
    const size_type extendedNumPressure = fracture->getNumExtendedPressure();

    for ( size_type i = 0; i < extendedNumPressure; ++i )
    {
        const size_type ii = extendedPressure [ i ];
        const base_node pointFlat = fracture->getMeshFEMPressure().point_of_basic_dof(ii);
        base_node pointMapped(0,0);
        pointMapped[0] = pointFlat[0];
        pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
        const scalar_type levelSetValue = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );

        if ( levelSetValue < 0 )
        {
            (*D) [ ii ] += DIn [ ii ];
            (*D) [ i + shiftPressure ] += DOut [ ii ];
        }
        else
        {
            (*D) [ ii ] += DOut [ ii ];
            (*D) [ i + shiftPressure ] += DIn [ ii ];
        }
    }

} // assembling_SourceF

void coupleFractures ( sparseMatrixPtr_Type& M,
                       const FracturesSetPtr_Type& fractures )
{
    const size_type numFractures = fractures->getNumberFractures ();

    for ( size_type i = 0; i < numFractures; ++i )
    {
        const FractureHandlerPtr_Type& fracture = fractures->getFracture(i);
        const pairSizeVectorContainer_Type& intersectElementsGlobalIndex =
                            fracture->getFractureIntersectElementsGlobalIndex ();

        for ( size_type j = i+1; j < numFractures; ++j )
        {
            const size_type numIntersections = intersectElementsGlobalIndex [j].size();
            for ( size_type k = 0; k < numIntersections; ++k )
            {
                const size_type first = intersectElementsGlobalIndex [j] [k].first;
                const size_type second = intersectElementsGlobalIndex [j] [k].second;
         
                (*M)( first, first ) = 1;
                (*M)( first, second ) = -1;
            }
        }
    }

} // coupleFractures

void velocityJump ( sparseMatrixPtr_Type& M,
                    const FractureHandlerPtr_Type& fracture,
                    const FractureHandlerPtr_Type& otherFracture,
                    const size_type& convex )
{
    const size_type shiftVelocity = fracture->getMeshFEMVelocity().nb_dof();
    scalarVector_Type V ( shiftVelocity, 0. );

    const size_type otherFractureId = otherFracture->getId();
    LevelSetHandlerPtr_Type& levelSetOtherFracture = otherFracture->getLevelSet();

    const getfem::pintegration_method intTypeIM = getfem::int_method_descriptor ( "IM_GAUSS1D(3)" );
    getfem::mesh_im_level_set meshImLevel ( *fracture->getMeshLevelSetIntersect ( otherFractureId ),
                                            getfem::mesh_im_level_set::INTEGRATE_BOUNDARY );

    meshImLevel.set_integration_method ( fracture->getMeshFlat().convex_index(), intTypeIM );

    meshImLevel.set_simplex_im ( intTypeIM );

    getfem::mesh_region meshElement;
    meshElement.add ( convex );

    generic_assembly assem;

    assem.set ( "V(#1)+=comp(vBase(#1).NonLin(#2))"
                "(:,1,1)" );

    level_set_unit_normal nterm ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem(),
                                  fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->values() );

    assem.push_mi ( meshImLevel );

    // Assign the mesh finite element space
    assem.push_mf ( fracture->getMeshFEMVelocity() );

    // Assign the mesh finite element space for the coefficients
    assem.push_mf ( fracture->getMeshLevelSetIntersect ( otherFractureId )->get_level_set(0)->get_mesh_fem() );

    // Assign the non linear term
    assem.push_nonlinear_term ( &nterm );

    // Set the matrices to save the evaluations
    assem.push_vec ( V );

    assem.assembly ( meshElement );

    const sizeVector_Type& extendedVelocity = fracture->getExtendedVelocity();
    const size_type extendedNumVelocity = fracture->getNumExtendedVelocity();

    for ( size_type i = 0; i < extendedNumVelocity; ++i )
    {
        const size_type ii = extendedVelocity [ i ];
        const base_node pointFlat = fracture->getMeshFEMVelocity().point_of_basic_dof(ii);
        base_node pointMapped(0,0);
        pointMapped[0] = pointFlat[0];
        pointMapped[1] = fracture->getLevelSet()->getData()->z_map( pointFlat );
        const scalar_type levelSetValue = levelSetOtherFracture->getData()->levelSetFunction ( pointMapped );

        if ( levelSetValue < 0 )
        {
            (*M) ( ii, 0 ) += V [ ii ];
            (*M) ( i + shiftVelocity, 0 ) -= V [ ii ];
        }
        else
        {
            (*M) ( ii, 0 ) -= V [ ii ];
            (*M) ( i + shiftVelocity, 0 ) += V [ ii ];
        }
    }

} // velocityJump

// L2 norm che tiene conto della frattura
scalar_type extended_L2_norm ( const scalarVectorPtr_Type& V1,
                               const LevelSetHandlerPtr_Type& levelSet,
                               const MeshHandlerPtr_Type& mesh,
                               const sizeVector_Type &eXt_dof,
                               const size_type& cutRegionFlag,
                               const size_type& uncutRegionFlag )
{

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    scalarVector_Type v(1), vIn(1), vOut(1);

    scalarVector_Type U(shiftVelocity, 0.), UIn(shiftVelocity, 0.);
    scalarVector_Type UOut(shiftVelocity, 0.);

    for ( size_type i = 0; i < eXt_dof.size(); ++i )
    {
        size_type ii = eXt_dof [ i ];
        if ( levelSet->getBaricenterValue(ii) < 0 )
        {
            UIn [ ii ] = (*V1) [ ii ];
            UOut [ ii ] = (*V1) [ i + shiftVelocity ];
        }
        else
        {
            UOut [ ii ] = (*V1) [ ii ];
            UIn [ ii ] = (*V1) [ i + shiftVelocity ];
        }
    }

    generic_assembly assem_Source, assem_SourceIn, assem_SourceOut;
    if ( mesh->getMeshFEMVector().get_qdim() == 1 )
    {
        assem_Source.set(
                "u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
        assem_SourceIn.set(
                "u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
        assem_SourceOut.set(
                "u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
    }
    else
    {
        assem_Source.set(
                "u=data(#1); V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k, j,k)");
        assem_SourceIn.set(
                "u=data(#1); V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k, j,k)");
        assem_SourceOut.set(
                "u=data(#1); V()+=u(i).u(j).comp(vBase(#1).vBase(#1))(i,k, j,k)");
    }

    assem_Source.push_mi(mesh->getIntegrationMethodVector());
    assem_SourceIn.push_mi(levelSet->getIntegrationMethodInside());
    assem_SourceOut.push_mi(levelSet->getIntegrationMethodOutside());

    assem_Source.push_mf(mesh->getMeshFEMVector());
    assem_SourceIn.push_mf(mesh->getMeshFEMVector());
    assem_SourceOut.push_mf(mesh->getMeshFEMVector());

    assem_Source.push_data(U);
    assem_SourceIn.push_data(UIn);
    assem_SourceOut.push_data(UOut);

    assem_Source.push_vec(v);
    assem_SourceIn.push_vec(vIn);
    assem_SourceOut.push_vec(vOut);

    assem_Source.assembly(mesh->getRegion(uncutRegionFlag));
    assem_SourceIn.assembly(mesh->getRegion(cutRegionFlag));
    assem_SourceOut.assembly(mesh->getRegion(cutRegionFlag));

    return std::sqrt(v [ 0 ] + vIn [ 0 ] + vOut [ 0 ]);
} // extended_L2_norm

// L2 norma della divergenza
scalar_type extended_divL2_norm ( const scalarVectorPtr_Type& V1,
                                  const LevelSetHandlerPtr_Type& levelSet,
                                  const MeshHandlerPtr_Type& mesh,
                                  const sizeVector_Type &eXt_dof,
                                  const size_type& cutRegionFlag,
                                  const size_type& uncutRegionFlag )
{

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    scalarVector_Type v(1), vIn(1), vOut(1);

    scalarVector_Type U(shiftVelocity), UIn(shiftVelocity, 0.), UOut(
            shiftVelocity, 0.);

    for ( size_type i = 0; i < eXt_dof.size(); ++i )
    {
        size_type ii = eXt_dof [ i ];
        if ( levelSet->getBaricenterValue(ii) < 0 )
        {
            UIn [ ii ] = (*V1) [ ii ];
            UOut [ ii ] = (*V1) [ i + shiftVelocity ];
        }
        else
        {
            UOut [ ii ] = (*V1) [ ii ];
            UIn [ ii ] = (*V1) [ i + shiftVelocity ];
        }
    }

    generic_assembly assem_Source, assem_SourceIn, assem_SourceOut;
    if ( mesh->getMeshFEMVector().get_qdim() == 1 )
    {
        assem_Source.set(
                "u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
        assem_SourceIn.set(
                "u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
        assem_SourceOut.set(
                "u=data(#1); V()+=u(i).u(j).comp(Base(#1).Base(#1))(i,j)");
    }
    else
    {
        assem_Source.set(
                "u=data(#1); V()+=u(i).u(j).comp(vGrad(#1).vGrad(#1))(i,k,k, j,l,l)");
        assem_SourceIn.set(
                "u=data(#1); V()+=u(i).u(j).comp(vGrad(#1).vGrad(#1))(i,k,k, j,l,l)");
        assem_SourceOut.set(
                "u=data(#1); V()+=u(i).u(j).comp(vGrad(#1).vGrad(#1))(i,k,k, j,l,l)");
    }

    assem_Source.push_mi(mesh->getIntegrationMethodVector());
    assem_SourceIn.push_mi(levelSet->getIntegrationMethodInside());
    assem_SourceOut.push_mi(levelSet->getIntegrationMethodOutside());

    assem_Source.push_mf(mesh->getMeshFEMVector());
    assem_SourceIn.push_mf(mesh->getMeshFEMVector());
    assem_SourceOut.push_mf(mesh->getMeshFEMVector());

    assem_Source.push_data(U);
    assem_SourceIn.push_data(UIn);
    assem_SourceOut.push_data(UOut);

    assem_Source.push_vec(v);
    assem_SourceIn.push_vec(vIn);
    assem_SourceOut.push_vec(vOut);

    assem_Source.assembly(mesh->getRegion(uncutRegionFlag));
    assem_SourceIn.assembly(mesh->getRegion(cutRegionFlag));
    assem_SourceOut.assembly(mesh->getRegion(cutRegionFlag));

    return sqrt(v [ 0 ] + vIn [ 0 ] + vOut [ 0 ]);
} // extended_divL2_norm

// calcola il flusso uscente da parte del bordo del dominio
void flux_out ( scalar_type& flux2d,
                const scalarVectorPtr_Type& V1,
                const LevelSetHandlerPtr_Type& levelSet,
                const MeshHandlerPtr_Type& mesh,
                const sizeVector_Type& eXt_dof,
                const sizeVector_Type& Boundary_cut,
                const sizeVector_Type& Boundary_uncut )
{

    const size_type shiftVelocity = mesh->getMeshFEMVector().nb_dof();
    scalarVector_Type v(1), vIn(1), vOut(1);

    scalarVector_Type U(shiftVelocity), UIn(shiftVelocity, 0.), UOut(
            shiftVelocity, 0.);

    for ( size_type i = 0; i < shiftVelocity; ++i )
    {
        U [ i ] = (*V1) [ i ];
    }

    for ( size_type i = 0; i < eXt_dof.size(); ++i )
    {
        size_type ii = eXt_dof [ i ];
        if ( levelSet->getDOFValue(ii) < 0 )
        {
            UIn [ ii ] = (*V1) [ ii ];
            UOut [ ii ] = (*V1) [ i + shiftVelocity ];
        }
        else
        {
            UOut [ ii ] = (*V1) [ ii ];
            UIn [ ii ] = (*V1) [ i + shiftVelocity ];
        }
    }

    generic_assembly assem_flux, assem_fluxIn, assem_fluxOut;
    if ( mesh->getMeshFEMVector().get_qdim() == 1 )
    {
        assem_flux.set("u=data$1(#1);  V()+=u(i).comp(Base(#1))(i)");
        assem_fluxIn.set("u=data$1(#1);  V()+=u(i).comp(Base(#1))(i)");
        assem_fluxOut.set("u=data$1(#1);  V()+=u(i).comp(Base(#1))(i)");
    }
    else
    {
        assem_flux.set("u=data(#1); V()+=u(i).comp(vBase(#1).Normal())(i,k,k)");
        assem_fluxIn.set(
                "u=data(#1); V()+=u(i).comp(vBase(#1).Normal())(i,k,k)");
        assem_fluxOut.set(
                "u=data(#1); V()+=u(i).comp(vBase(#1).Normal())(i,k,k)");
    }

    assem_flux.push_mi(mesh->getIntegrationMethodVector());
    assem_fluxIn.push_mi(levelSet->getIntegrationMethodInside());
    assem_fluxOut.push_mi(levelSet->getIntegrationMethodOutside());

    assem_flux.push_mf(mesh->getMeshFEMVector());
    assem_fluxIn.push_mf(mesh->getMeshFEMVector());
    assem_fluxOut.push_mf(mesh->getMeshFEMVector());

    assem_flux.push_data(U);
    assem_fluxIn.push_data(UIn);
    assem_fluxOut.push_data(UOut);

    assem_flux.push_vec(v);
    assem_fluxIn.push_vec(vIn);
    assem_fluxOut.push_vec(vOut);

    for ( size_type bndID = 0; bndID < Boundary_uncut.size(); bndID++ )
    {
        assem_flux.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(
                        Boundary_uncut [ bndID ]));
    }
    for ( size_type bndID = 0; bndID < Boundary_cut.size(); bndID++ )
    {
        assem_fluxIn.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(
                        Boundary_cut [ bndID ]));
        assem_fluxOut.assembly(
                mesh->getMeshFEMVector().linked_mesh().get_mpi_sub_region(
                        Boundary_cut [ bndID ]));
    }

    flux2d += v [ 0 ] + vIn [ 0 ] + vOut [ 0 ];
} // flux_out

// L2 norm con peso
scalar_type extended_L2_norm ( const scalarVectorPtr_Type& V1,
                               const scalarVector_Type& weight,
                               const LevelSetHandlerPtr_Type& levelSet,
                               const mesh_im& mediumIntegrationMethodVelocity,
                               const mesh_fem& mf,
                               const mesh_fem& mf2,
                               const sizeVector_Type& eXt_dof,
                               const mesh_region& cutRegion,
                               const mesh_region& uncutRegion )
{
    size_type dof_shift = mf.nb_dof();
    scalarVector_Type v(mf2.nb_dof(), 0.), vIn(mf2.nb_dof(), 0.), vOut(
            mf2.nb_dof(), 0.);
    scalar_type vv, vvIn, vvOut;

    scalarVector_Type U(mf.nb_dof()), UIn(mf.nb_dof(), 0.), UOut(mf.nb_dof(),
            0.);

    for ( size_type i = 0; i < eXt_dof.size(); ++i )
    {
        size_type ii = eXt_dof [ i ];
        if ( levelSet->getBaricenterValue(ii) < 0 )
        {
            UIn [ ii ] = (*V1) [ ii ];
            UOut [ ii ] = (*V1) [ i + dof_shift ];
        }
        else
        {
            UOut [ ii ] = (*V1) [ ii ];
            UIn [ ii ] = (*V1) [ i + dof_shift ];
        }
    }

    generic_assembly assem_Source, assem_SourceIn, assem_SourceOut;
    if ( mf.get_qdim() == 1 )
    {
        assem_Source.set("u=data(#2);"
            "V(#1)+=u(i).u(j).comp(Base(#2).Base(#2).Base(#1))(i,j,:)");
        assem_SourceIn.set("u=data(#2);"
            "V(#1)+=u(i).u(j).comp(Base(#2).Base(#2).Base(#1))(i,j,:)");
        assem_SourceOut.set("u=data(#2);"
            "V(#1)+=u(i).u(j).comp(Base(#2).Base(#2).Base(#1))(i,j,:)");
    }
    else
    {
        assem_Source.set(
                "u=data(#2); V(#1)+=u(i).u(j).comp(vBase(#2).vBase(#2).Base(#1))(i,k, j,k,:)");
        assem_SourceIn.set(
                "u=data(#2); V(#1)+=u(i).u(j).comp(vBase(#2).vBase(#2).Base(#1))(i,k, j,k,:)");
        assem_SourceOut.set(
                "u=data(#2); V(#1)+=u(i).u(j).comp(vBase(#2).vBase(#2).Base(#1))(i,k, j,k,:)");
    }

    assem_Source.push_mi(mediumIntegrationMethodVelocity);
    assem_SourceIn.push_mi(levelSet->getIntegrationMethodInside());
    assem_SourceOut.push_mi(levelSet->getIntegrationMethodOutside());

    assem_Source.push_mf(mf2);
    assem_SourceIn.push_mf(mf2);
    assem_SourceOut.push_mf(mf2);

    assem_Source.push_mf(mf);
    assem_SourceIn.push_mf(mf);
    assem_SourceOut.push_mf(mf);

    assem_Source.push_data(U);
    assem_SourceIn.push_data(UIn);
    assem_SourceOut.push_data(UOut);

    assem_Source.push_vec(v);
    assem_SourceIn.push_vec(vIn);
    assem_SourceOut.push_vec(vOut);

    assem_Source.assembly(uncutRegion);
    assem_SourceIn.assembly(cutRegion);
    assem_SourceOut.assembly(cutRegion);

    vv = 0;
    vvIn = 0;
    vvOut = 0;
    for ( size_type i = 0; i < mf2.nb_dof(); ++i )
    {
        vv += v [ i ] * weight [ i ];
        vvIn += vIn [ i ] * weight [ i ];
        vvOut += vOut [ i ] * weight [ i ];
    }

    return sqrt(vv + vvIn + vvOut);

} // extended_L2_norm

}// namespace getfem
