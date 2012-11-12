package src;

public class PopulatorDerivator<T> implements IDerivator<T>
{
	private final IDerivator<T> derivator;
	private final IPopulator<? super T> populator;
	
	public PopulatorDerivator(IDerivator<T> derivator, IPopulator<? super T> populator)
	{
		this.derivator = derivator;
		this.populator = populator;
	}
	
	public T derive(T values)
	{
		T result = derivator.derive(values);
		populator.populate(result);
		return result;
	}
}
