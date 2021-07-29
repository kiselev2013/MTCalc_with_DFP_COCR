#pragma once
#include "PointVector.h"

/*! типы результатов, которые могут быть выданы 3d сетками */
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

/*! базовый класс дл€ 3d  Ё сеток (с геометрически однотипными элементами) */
class AbstractFEM3D
{
public:
	/*! кол-во узлов */
	virtual int GetNumberOfNodes()=0;

	/*! кол-во элементов */
	virtual int GetNumberOfElements()=0;

	/*! кол-во узлов элемента */
	virtual int GetElementNodesNumber()=0;

	/*! координаты узла по номеру */
	virtual const pv::Point3D GetNode(const int& i_node)=0;

	/*! координаты узла по номеру */
	virtual const pv::Point3D GetNodeTrue(const int& i_node)=0;

	/*! номер узла в элементе */
	virtual int GetNodeNumberOnElement(const int& i_element, const int& i_node)=0;

	/*! номер материала элемента */
	virtual int GetElementMaterial(const int& i_element)=0;

	/*! тип элемента */
	virtual int GetTypeOfElement(const int& i_element)=0;

	/*! значение результата в центре элемента */
	virtual double GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type)=0;

	/*! кол-во точек выдачи результата */
	virtual int GetNumberOfResPoints(const Res3DValueType& r_type)=0;

	/*! точка выдачи результата */
	virtual pv::Point3D GetResPoint(const Res3DValueType& r_type, const int& i_point)=0;

	/*! указатель на массив regular */
	virtual int * GetPointerToRegular()=0;

	/*! размерность массива X регул€рной сетки */
	virtual int GetXSize()=0;

	/*! размерность массива Y регул€рной сетки */
	virtual int GetYSize()=0;

	/*! размерность массива Z регул€рной сетки */
	virtual int GetZSize()=0;

	/*! массив X регул€рной сетки */
	virtual double *GetPointerToX()=0;

	/*! массив Y регул€рной сетки */
	virtual double *GetPointerToY()=0;

	/*! массив Z регул€рной сетки */
	virtual double *GetPointerToZ()=0;

	/*! сохранить результат */
	virtual void SaveResult(const Res3DValueType& r_type, const double& r_value, const int& i_point, const int& i_time)=0;
};

