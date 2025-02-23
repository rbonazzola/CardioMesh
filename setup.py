from setuptools import setup, find_packages

setup(
    name="cardio-mesh",
    version="0.0.1",
    author="Rodrigo Bonazzola (rbonazzola)",
    author_email="rodbonazzola@gmail.com",
    description="Python package for handling cardiac meshes",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/rbonazzola/CardioMesh",
    packages=["cardio_mesh"],
    install_requires=[
        "numpy>=1.24.0",
        "meshio>=5.3.5",
        "vtk>=9.3.1",
        "vedo>2024.5.2",
        "scipy>=1.10.1",
        "trimesh>=4.5.3"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
