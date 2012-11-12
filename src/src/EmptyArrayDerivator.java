package src;

public class EmptyArrayDerivator implements IDerivator<Array2D>
{
	@Override
	public Array2D derive(Array2D values)
	{
		return new Array2D(values.width, values.height);
	}

}
