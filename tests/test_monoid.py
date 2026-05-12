import pytest
import sys
import os

# Добавляем корневую папку проекта в путь (если не установлен пакет)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src._core.monoid import Monoid
# from abc import ABC

class ConcreteMonoid(Monoid):
    """Конкретный наследник для тестирования (не абстрактный)."""
    def __init__(self, components, dual=False, multiplicity=1):
        super().__init__(components, dual, multiplicity)
    def __mul__(self, other):
        # Минимальная реализация для тестирования базового поведения
        result = super().__mul__(other)
        if result is not NotImplemented: return result
        # Имитация умножения двух конкретных моноидов (просто возвращаем себя)
        return self

def test_zero():
    z = ConcreteMonoid.zero()
    assert z.is_zero()
    assert not z.is_one()
    assert len(z.components) == 0
    assert z.multiplicity == -1

def test_one():
    o = ConcreteMonoid.one()
    assert o.is_one()
    assert not o.is_zero()
    assert len(o.components) == 0
    assert o.multiplicity == -2

def test_zero_one_equality():
    z = ConcreteMonoid.zero()
    o = ConcreteMonoid.one()
    assert z != o
    assert z == ConcreteMonoid.zero()
    assert o == ConcreteMonoid.one()

def test_hash():
    z1 = ConcreteMonoid.zero()
    z2 = ConcreteMonoid.zero()
    o1 = ConcreteMonoid.one()
    o2 = ConcreteMonoid.one()
    assert hash(z1) == hash(z2)
    assert hash(o1) == hash(o2)
    assert hash(z1) != hash(o1)

def test_multiplication_with_zero():
    z = ConcreteMonoid.zero()
    any_m = ConcreteMonoid([1,2])
    assert (z * any_m).is_zero()
    assert (any_m * z).is_zero()
    assert (z * z).is_zero()

def test_multiplication_with_one():
    o = ConcreteMonoid.one()
    m = ConcreteMonoid([1,2])
    assert (o * m) is m
    assert (m * o) is m