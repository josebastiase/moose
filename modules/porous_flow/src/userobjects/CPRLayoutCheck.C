//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CPRLayoutCheck.h"
#include "FEProblemBase.h"
#include "NonlinearSystemBase.h"
#include "MooseMesh.h"

#include "libmesh/system.h"
#include "libmesh/dof_map.h"
#include "libmesh/dof_object.h"
#include "libmesh/fe_type.h"

registerMooseObject("PorousFlowApp", CPRLayoutCheck);

InputParameters
CPRLayoutCheck::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription(
      "Verifies that the dofs of all the nonlinear variables are interleaved per node (finite "
      "elements) or per element (finite volumes), as required by '-pc_type cpr', and prints the "
      "PETSc options to use.  Run with --node-major-dofs on the command line to get this layout");
  params.addRequiredParam<NonlinearVariableName>(
      "pressure_variable", "Nonlinear variable treated as the pressure by CPR (eg porepressure)");
  params.set<ExecFlagEnum>("execute_on") = EXEC_INITIAL;
  params.suppressParameter<ExecFlagEnum>("execute_on");
  return params;
}

CPRLayoutCheck::CPRLayoutCheck(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _pressure_var(getParam<NonlinearVariableName>("pressure_variable"))
{
}

void
CPRLayoutCheck::initialSetup()
{
  const libMesh::System & sys = _fe_problem.getNonlinearSystemBase(/*nl_sys_num=*/0).system();
  const unsigned int sn = sys.number();
  const unsigned int nvars = sys.n_vars();

  if (!sys.has_variable(_pressure_var))
    paramError("pressure_variable", "'", _pressure_var, "' is not a nonlinear variable");
  const unsigned int pfield = sys.variable_number(_pressure_var);

  // 0 = ok, 1 = a variable does not have one dof per node/element, 2 = not interleaved,
  // 3 = dofs that belong to no node or element (eg scalar variables)
  unsigned int bad = 0;
  dof_id_type counted_dofs = 0;
  std::string where;

  // Every node (finite elements) or element (finite volumes) that holds dofs of this system must
  // hold one dof of every variable, numbered contiguously in the order of the variables
  const auto check = [&](const libMesh::DofObject & obj, const std::string & what)
  {
    if (obj.n_dofs(sn) == 0)
      return;
    if (obj.n_dofs(sn) != nvars)
    {
      bad = std::max(bad, 1u);
      where = what;
      return;
    }
    const dof_id_type d0 = obj.dof_number(sn, 0, 0);
    for (const auto v : make_range(nvars))
    {
      if (obj.n_comp(sn, v) != 1)
      {
        bad = std::max(bad, 1u);
        where = what;
      }
      else if (obj.dof_number(sn, v, 0) != d0 + v)
      {
        bad = std::max(bad, 2u);
        where = what;
      }
    }
    counted_dofs += nvars;
  };

  // Scalar variables have dofs that belong to no node or element (libMesh stores them on the last
  // node, which would otherwise be reported as a node with the wrong number of dofs)
  for (const auto v : make_range(nvars))
    if (sys.variable_type(v).family == libMesh::SCALAR)
      bad = 3;

  const auto & mesh = _fe_problem.mesh().getMesh();
  if (bad == 0)
  {
    for (const auto & node : mesh.local_node_ptr_range())
      check(*node, "node " + std::to_string(node->id()));
    for (const auto & elem : mesh.active_local_element_ptr_range())
      check(*elem, "element " + std::to_string(elem->id()));
    if (bad == 0 && counted_dofs != sys.get_dof_map().n_local_dofs())
      bad = 3;
  }

  _communicator.max(bad);

  switch (bad)
  {
    case 1:
      mooseError(name(),
                 ": every nonlinear variable must have exactly one dof on every node (finite "
                 "elements) or element (finite volumes) where the nonlinear variables live, "
                 "eg first-order LAGRANGE or finite-volume variables on the same blocks.  Problem "
                 "at ",
                 where.empty() ? "some node or element" : where,
                 ".  Use PCCPRSetNodeDofs() for general layouts");
    case 2:
      mooseError(name(),
                 ": the dofs are not interleaved per node or element (",
                 where,
                 ").  Run with --node-major-dofs on the command line, and make sure all the "
                 "nonlinear variables have the same type and blocks");
    case 3:
      mooseError(name(),
                 ": the nonlinear system has dofs that belong to no node or element (eg scalar "
                 "variables), so block-size CPR cannot be used");
    default:
      break;
  }

  _console << "\n"
           << name() << ": dof layout OK for CPR.  Use\n"
           << "  -pc_type cpr -pc_cpr_block_size " << nvars << " -pc_cpr_pressure_field " << pfield
           << "\n  (pressure variable '" << _pressure_var << "' is field " << pfield << " of "
           << nvars << ")\n"
           << std::endl;
}
