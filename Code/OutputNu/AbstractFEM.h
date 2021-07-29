#pragma once
#include "PointVector.h"

/*! ���� �����������, ������� ����� ���� ������ 3d ������� */
enum Res3DValueType
{
	vtWithDiscontinuity,
	vtWithoutDiscontinuity,
	vtU,
	vtdUdx,
	vtdUdy,
	vtdUdz,
	vtExSin,
	vtExCos,
	vtEySin,
	vtEyCos,
	vtEzSin,
	vtEzCos,
	vtAx,
	vtAy,
	vtAz,
	vtAxSin,
	vtAxCos,
	vtAySin,
	vtAyCos,
	vtAzSin,
	vtAzCos,
	vtRotxA,
	vtRotyA,
	vtRotzA,
	vtRotxASin,
	vtRotxACos,
	vtRotyASin,
	vtRotyACos,
	vtRotzASin,
	vtRotzACos
};

/*! ������� ����� ��� 3d �� ����� (� ������������� ����������� ����������) */
class AbstractFEM3D
{
public:
	/*! ���-�� ����� */
	virtual int GetNumberOfNodes()=0;

	/*! ���-�� ��������� */
	virtual int GetNumberOfElements()=0;

	/*! ���-�� ����� �������� */
	virtual int GetElementNodesNumber()=0;

	/*! ���������� ���� �� ������ */
	virtual const pv::Point3D GetNode(const int& i_node)=0;

	/*! ���������� ���� �� ������ */
	virtual const pv::Point3D GetNodeTrue(const int& i_node)=0;

	/*! ����� ���� � �������� */
	virtual int GetNodeNumberOnElement(const int& i_element, const int& i_node)=0;

	/*! ����� ��������� �������� */
	virtual int GetElementMaterial(const int& i_element)=0;

	/*! ��� �������� */
	virtual int GetTypeOfElement(const int& i_element)=0;

	/*! �������� ���������� � ������ �������� */
	virtual double GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type)=0;

	/*! ���-�� ����� ������ ���������� */
	virtual int GetNumberOfResPoints(const Res3DValueType& r_type)=0;

	/*! ����� ������ ���������� */
	virtual pv::Point3D GetResPoint(const Res3DValueType& r_type, const int& i_point)=0;

	/*! ��������� �� ������ regular */
	virtual int * GetPointerToRegular()=0;

	/*! ����������� ������� X ���������� ����� */
	virtual int GetXSize()=0;

	/*! ����������� ������� Y ���������� ����� */
	virtual int GetYSize()=0;

	/*! ����������� ������� Z ���������� ����� */
	virtual int GetZSize()=0;

	/*! ������ X ���������� ����� */
	virtual double *GetPointerToX()=0;

	/*! ������ Y ���������� ����� */
	virtual double *GetPointerToY()=0;

	/*! ������ Z ���������� ����� */
	virtual double *GetPointerToZ()=0;

	/*! ��������� ��������� */
	virtual void SaveResult(const Res3DValueType& r_type, const double& r_value, const int& i_point, const int& i_time)=0;
};

